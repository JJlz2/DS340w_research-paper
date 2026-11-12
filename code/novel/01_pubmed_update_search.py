"""Novel approach, step 1 of 4 (Acquisition).

Reruns the parent paper's own PubMed search string (Amdi & King, 2025,
Supplementary file 1) for records added to PubMed after the parent's
search date, through the NCBI E-utilities API. The output is the list of
candidate studies a living update of the meta-analysis has to screen.

The parent searched "from inception to 10 November 2023" and kept search
alerts running "until the week before submission" (submitted ~23 September
2024), so records added after mid-September 2024 were never seen by it.

Run from the repo root:
    python code/novel/01_pubmed_update_search.py
Output: data/processed/pubmed_update_search.csv
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
# The parent's PubMed string, verbatim except straight quotes
QUERY = ('(eccentric) AND (duration OR phase OR speed OR cadence OR temp*) '
         'AND ("resistance training" OR "strength training")')
PARENT_SEARCH_END = "2023/11/10"   # main search end date
ALERTS_END = "2024/09/16"          # approximate end of the parent's search alerts
# Studies screened in this project so far (PMID -> decision)
SCREENED = {"42507730": "eligible - Hao et al. (2026), added to the meta-analysis"}

OUT = Path("data/processed/pubmed_update_search.csv")


def get_json(endpoint, params):
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        return json.load(resp)


def main():
    search = get_json("esearch.fcgi", {
        "db": "pubmed", "term": QUERY, "retmode": "json", "retmax": 10000,
        "datetype": "edat", "mindate": PARENT_SEARCH_END, "maxdate": "3000"})
    pmids = search["esearchresult"]["idlist"]

    rows = []
    for start in range(0, len(pmids), 200):
        batch = pmids[start:start + 200]
        summary = get_json("esummary.fcgi", {
            "db": "pubmed", "id": ",".join(batch), "retmode": "json"})["result"]
        for pmid in batch:
            rec = summary[pmid]
            entrez = next((h["date"][:10] for h in rec.get("history", [])
                           if h.get("pubstatus") == "entrez"), "")
            rows.append({
                "pmid": pmid,
                "first_author": rec.get("sortfirstauthor", ""),
                "pub_date": rec.get("pubdate", ""),
                "title": rec.get("title", ""),
                "journal": rec.get("fulljournalname", ""),
                "entrez_date": entrez,
                "after_parent_alerts": entrez.replace("-", "/") > ALERTS_END,
                "screening_status": SCREENED.get(pmid, "not yet screened"),
            })
        time.sleep(0.4)   # NCBI asks for <= 3 requests per second without a key

    rows.sort(key=lambda r: r["entrez_date"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    after = sum(r["after_parent_alerts"] for r in rows)
    print(f"{len(rows)} records added to PubMed since {PARENT_SEARCH_END}; "
          f"{after} after the parent's search alerts ended (~{ALERTS_END}).")
    print(f"Screened so far: {len(SCREENED)}. Wrote {OUT}")


if __name__ == "__main__":
    main()
