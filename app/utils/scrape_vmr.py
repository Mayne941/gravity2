import os
import re
import requests
import struct
import pandas as pd
from bs4 import BeautifulSoup as bs


class Scraper:
    '''Get latest VMR from ICTV, structure to work with GRAViTy'''

    def __init__(self, dir) -> None:
        self.url = "https://ictv.global/vmr"
        self.url_stem = "https://ictv.global"
        self.save_dir = dir
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
            # RM < Check the following are correct
            "ssDNA(+)": "II",
            "ssDNA(-)": "II",
            "ssDNA(+/-)": "II",
            "ssRNA": "???",
            "ssRNA(+/-)": "???",

        }

    def scrape(self) -> pd.DataFrame:
        '''Scrape URL, extract links table, get latest link, parse as excel table'''
        soup = bs(requests.get(self.url).text, "html.parser")
        table_tags = soup.find("table", {"id": "edit-table"})
        url_rows = [re.search(r"<a\shref=.*\">", str(row))[0].replace('<a href=\"', '').replace(
            '">', '') for row in table_tags.findAll("tr") if re.search(r"<a\shref=.*</a>", str(row))]
        latest_url = f"{self.url_stem}{url_rows[0]}"
        return pd.read_excel(latest_url)

    def get_baltimore(self, row):
        '''Make new column for roman numeral notation Baltimore classification.'''
        try:
            return self.baltimore_map[row["Genome composition"]]
        except KeyError:
            print(
                f"VMR row with unknown genome composition detected ({row['Genome composition']}). Skipping row...")
            return ""

    def etl(self, df) -> pd.DataFrame:
        '''Extract, transform, load.'''
        df["Baltimore Group"] = df.apply(
            lambda x: self.get_baltimore(x), axis=1)
        df["Genetic code table"] = ""  # RM < What should this be?
        df["Taxonomic grouping"] = ""  # RM < What should this be?
        return df

    def save_csv(self, df) -> None:
        df.to_csv(f"{self.save_dir}latest_vmr.csv")

    def main(self) -> str:
        '''Entrypoint'''
        links_table = self.scrape()
        structured_table = self.etl(links_table)
        self.save_csv(structured_table)
        return self.save_dir


def scrape(payload):
    s = Scraper(payload["save_path"])
    try:
        dir = s.main()
        return f"Success! VMR saved to {dir}"
    except:
        print("Whoops")
        return f"Scrape unsuccessful: possibly the VMR format and/or URL has changed.\n Please check url and table formatting in app/utils/scrape_vmr.py, or contact your administrator."


if __name__ == "__main__":
    scrape("./data/")
