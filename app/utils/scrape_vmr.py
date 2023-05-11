import os
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup as bs

from app.utils.timer import timing


class Scraper:
    '''Get latest VMR from ICTV, structure to work with GRAViTy'''

    def __init__(self, payload) -> None:
        self.url = "https://ictv.global/vmr"
        self.url_stem = "https://ictv.global"
        self.save_dir = payload["save_path"]
        self.fname = payload["vmr_name"]
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)
        self.baltimore_map = {
            "dsDNA": "I",
            "ssDNA": "II",
            "dsRNA": "III",
            "ssRNA(+)": "IV",
            "ssRNA(-)": "V",
            "ssRNA-RT": "VI",
            "dsDNA-RT": "VII",
            "ssDNA(+)": "II",
            "ssDNA(-)": "II",
            "ssDNA(+/-)": "II",
            "ssRNA": "IV, V",
            "ssRNA(+/-)": "IV, V"
        }

    def scrape(self) -> pd.DataFrame:
        '''Scrape URL, extract links table, get latest link, parse as excel table'''
        soup = bs(requests.get(self.url).text, "html.parser")
        table_tags = soup.find("table", {"id": "edit-table"})
        url_rows = [re.search(r"<a\shref=.*\">", str(row))[0].replace('<a href=\"', '').replace(
            '">', '') for row in table_tags.findAll("tr") if re.search(r"<a\shref=.*</a>", str(row))]
        latest_url = f"{self.url_stem}{url_rows[1]}"
        return pd.read_excel(latest_url)

    def get_baltimore_and_code_table(self, row) -> pd.Series:
        '''Make new column for roman numeral notation Baltimore classification and genome coverage indicator.'''
        try:
            baltimore = self.baltimore_map[row["Genome composition"]]
        except KeyError:
            print(
                f"VMR row with unknown genome composition detected ({row['Genome composition']}). Skipping row...")
            baltimore = ""

        if "complete" in str(row["Genome coverage"]).lower():
            code_table = 1
        else:
            code_table = 0

        isolate_col = str(row["Virus isolate designation"]
                          ).replace("\n", "").replace(" ", "")

        return pd.Series([baltimore, code_table, isolate_col])

    def is_segmented(self, row) -> int:
        if ":" in row["Virus GENBANK accession"] or ";" in row["Virus GENBANK accession"]:
            return 1
        else:
            return 0

    def etl(self, df) -> pd.DataFrame:
        '''Extract, transform, load. Get baltimore group and synthesise column for genome coverage'''
        df[["Baltimore Group", "Genetic code table", "Virus isolate designation"]] = df.apply(
            lambda x: self.get_baltimore_and_code_table(x), axis=1)
        df["Taxonomic grouping"] = ""
        df = df.replace(",", ";", regex=True)

        '''Find, print list of and remove duplicates or entries with no Acc ID'''
        df[df["Virus GENBANK accession"].isin(df["Virus GENBANK accession"][df["Virus GENBANK accession"].duplicated(
        )])]["Virus GENBANK accession"].to_csv(f"{self.save_dir}/bad_accession.csv")
        df = df.drop_duplicates(subset="Virus GENBANK accession", keep=False)

        '''Exclude partial & segmented genomes'''
        df = df[df["Genetic code table"] == 1]

        return df

    def save_csv(self, df) -> None:
        df.to_csv(f"{self.save_dir}/{self.fname}")

    def main(self) -> str:
        '''Entrypoint'''
        links_table = self.scrape()
        structured_table = self.etl(links_table)
        self.save_csv(structured_table)
        return self.save_dir


def construct_first_pass_set(row, df, threshold) -> pd.Series:
    '''Construct first-pass grouping; ensure diversity is maintained in representative examples'''
    if df[df["Family"] == row["Family"]].shape[0] >= threshold:
        return row["Family"]
    else:
        if df[df["Genus"] == row["Genus"]].shape[0] >= threshold:
            return row["Genus"]
        else:
            return row["Species"]


@timing
def scrape(payload) -> str:
    print("VMR scrape started. This may take a short while.")
    s = Scraper(payload)
    try:
        dir = s.main()
        return f"Success! VMR saved to {dir}"
    except Exception as ex:
        print(f"Whoops: {ex}")
        return f"Scrape unsuccessful: possibly the VMR format and/or URL has changed.\n Please check url and table formatting in app/utils/scrape_vmr.py, or contact your administrator."


@timing
def first_pass_baltimore_filter(payload) -> str:
    try:
        df = pd.read_csv(f"{payload['save_path']}/{payload['vmr_name']}")
        df["Taxonomic grouping"] = df.apply(
            lambda x: construct_first_pass_set(x, df, payload["filter_threshold"]), axis=1)
        df = df.drop_duplicates(subset="Taxonomic grouping", keep="first")

        if payload["baltimore_filter"] == "RNA":
            df = df[(df["Baltimore Group"].str.contains("III")) | (df["Baltimore Group"].str.contains("IV")) | (
                df["Baltimore Group"].str.contains("V")) | (df["Baltimore Group"].str.contains("VI"))]
        elif payload["baltimore_filter"] == "dsRNA":
            df = df[(df["Baltimore Group"].str.contains("I")) |
                    (df["Baltimore Group"].str.contains("VII"))]
        elif payload["baltimore_filter"] == "ssDNA":
            df = df[df["Baltimore Group"].str.contains("II")]

        df.to_csv(f"{payload['save_path']}/{payload['save_name']}")
        return f"Success! VMR saved to ./{payload['save_path']}"
    except Exception as e:
        print(f"Error processing VMR, error: {e}")
        return f"Unsuccessful.\n Please check url and table formatting in app/utils/scrape_vmr.py, your spelling of the filter parameters, or contact your administrator."

@timing
def first_pass_taxon_filter(payload) -> str:
    try:
        df = pd.read_csv(f"{payload['save_path']}/{payload['vmr_name']}")
        df = df[df[payload['filter_level']] == payload["filter_name"]]
        df["Taxonomic grouping"] = df.apply(
            lambda x: construct_first_pass_set(x, df, payload["filter_threshold"]), axis=1)
        df = df.drop_duplicates(subset="Taxonomic grouping", keep="first")
        df.to_csv(f"{payload['save_path']}/{payload['save_name']}")
        return f"Success! VMR saved to ./{payload['save_path']}"
    except Exception as e:
        print(f"Error processing VMR, error: {e}")
        return f"Unsuccessful.\n Please check url and table formatting in app/utils/scrape_vmr.py, your spelling of the filter parameters, or contact your administrator."


@timing
def second_pass(payload) -> str:
    try:
        df = pd.read_csv(f"{payload['save_path']}/{payload['vmr_name']}")
        df = df[df[payload['filter_level']] == payload["filter_name"]]
        df.to_csv(f"{payload['save_path']}/{payload['save_name']}")
        return f"Success! VMR saved to ./{payload['save_path']}"
    except Exception as e:
        print(f"Error processing VMR, error: {e}")
        return f"Unsuccessful.\n Please check url and table formatting in app/utils/scrape_vmr.py, or contact your administrator."

@timing
def vmr_filter(payload) -> str:
    try:
        df = pd.read_csv(f"{payload['save_path']}/{payload['vmr_name']}")
        df = df.drop_duplicates(subset=payload["filter_level"], keep="first")
        df.to_csv(f"{payload['save_path']}/{payload['save_name']}")
        return f"Success! VMR saved to ./{payload['save_path']}"
    except Exception as e:
        print(f"Error processing VMR, error: {e}")
        return f"Unsuccessful.\n Please check url and table formatting in app/utils/scrape_vmr.py, or contact your administrator."


if __name__ == "__main__":
    scrape("./data/")
