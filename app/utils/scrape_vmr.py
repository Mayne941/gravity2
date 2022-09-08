import os
import re
import requests
import struct
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
        self.threshold = payload["filter_threshold"]
        self.first_pass_filter = payload["first_pass_filter"]
        self.second_pass_filter = payload["second_pass_filter"]
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
            "ssRNA": "II, III",
            "ssRNA(+/-)": "II, III"
        }

    def scrape(self) -> pd.DataFrame:
        '''Scrape URL, extract links table, get latest link, parse as excel table'''
        soup = bs(requests.get(self.url).text, "html.parser")
        table_tags = soup.find("table", {"id": "edit-table"})
        url_rows = [re.search(r"<a\shref=.*\">", str(row))[0].replace('<a href=\"', '').replace(
            '">', '') for row in table_tags.findAll("tr") if re.search(r"<a\shref=.*</a>", str(row))]
        latest_url = f"{self.url_stem}{url_rows[0]}"
        return pd.read_excel(latest_url)

    def get_baltimore_and_code_table(self, row):
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

    def construct_first_pass_set(self, row, df):
        '''Construct first-pass grouping; ensure diversity is maintained in representative examples'''
        # RM < Vectorise this as it's slow
        if df[df["Family"] == row["Family"]].shape[0] >= self.threshold:
            return row["Family"]
        else:
            if df[df["Genus"] == row["Genus"]].shape[0] >= self.threshold:
                return row["Genus"]
            else:
                return row["Species"]

    def etl(self, df) -> pd.DataFrame:
        '''Extract, transform, load. Get baltimore group and synthesise column for genome coverage'''
        df[["Baltimore Group", "Genetic code table", "Virus isolate designation"]] = df.apply(
            lambda x: self.get_baltimore_and_code_table(x), axis=1)
        # Not required for PL2 so default to blank
        df["Taxonomic grouping"] = ""

        '''Find, print list of and remove duplicates or entries with no Acc ID'''
        df[df["Virus GENBANK accession"].isin(df["Virus GENBANK accession"][df["Virus GENBANK accession"].duplicated(
        )])]["Virus GENBANK accession"].to_csv(f"{self.save_dir}/bad_accession.csv")
        df = df.drop_duplicates(subset="Virus GENBANK accession", keep=False)

        '''Exclude partial genomes'''
        df = df[df["Genetic code table"] == 1]

        if self.first_pass_filter:
            df["Taxonomic grouping"] = df.apply(
                lambda x: self.construct_first_pass_set(x, df), axis=1)
            df = df.drop_duplicates(subset="Taxonomic grouping", keep="first")
            print(
                f"Value counts for first pass filter, by Baltimore: {df['Baltimore Group'].value_counts()}")

        if self.second_pass_filter:
            ...

        return df

    def save_csv(self, df) -> None:
        df.to_csv(f"{self.save_dir}/{self.fname}")

    def main(self) -> str:
        '''Entrypoint'''
        links_table = self.scrape()
        structured_table = self.etl(links_table)
        self.save_csv(structured_table)
        return self.save_dir


@timing
def scrape(payload):
    print("VMR scrape started. This may take a short while.")
    s = Scraper(payload)
    try:
        dir = s.main()
        return f"Success! VMR saved to {dir}"
    except Exception as ex:
        print(f"Whoops: {ex}")
        return f"Scrape unsuccessful: possibly the VMR format and/or URL has changed.\n Please check url and table formatting in app/utils/scrape_vmr.py, or contact your administrator."


if __name__ == "__main__":
    scrape("./data/")
