"""Novel approach, step 1 of 4 (Acquisition).

Reruns the parent paper's own PubMed search string (Amdi & King, 2025,
Supplementary file 1) for records added to PubMed after the parent's
search date, through the NCBI E-utilities API. The output is the list of
candidate studies a living update of the meta-analysis has to screen.

The parent searched "from inception to 10 November 2023" and kept search
alerts running "until the week before submission" (submitted ~23 September
2024), so records added after mid-September 2024 were never seen by it.

The search used in the paper was run on 12 November 2026 and is saved as
data/processed/pubmed_update_search.csv. A rerun searches the same dates,
writes a separate file so the saved search is never overwritten, and lists
any records that differ from it. PubMed adds subject terms (MeSH) to records
weeks after they first appear, so an older record can start matching the
search later, and a rerun can find a few more records than the saved one.

Run from the repo root (needs internet; Python 3, standard library only):
    python code/novel/01_pubmed_update_search.py
Output: data/processed/pubmed_update_search_rerun.csv
"""
import csv
import json
import sys
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
SEARCH_END = "2026/11/12"          # the day the saved search (the one in the paper) was run
# Studies screened in this project so far (PMID -> decision)
SCREENED = {"42507730": "eligible - Hao et al. (2026), added to the meta-analysis"}

SAVED = Path("data/processed/pubmed_update_search.csv")       # the search used in the paper
OUT = Path("data/processed/pubmed_update_search_rerun.csv")   # a rerun never overwrites it
# Plain-English column names for the saved file
PLAIN = {"pmid": "PMID", "first_author": "First author", "pub_date": "Publication date",
         "title": "Title", "journal": "Journal", "entrez_date": "Date added to PubMed",
         "after_parent_alerts": "Added after the parent's search alerts ended",
         "screening_status": "Screening status"}


def get_json(endpoint, params, tries=4):
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except OSError as err:   # NCBI sometimes returns a one-off server error; wait and retry
            if attempt == tries:
                raise
            print(f"PubMed did not answer ({err}); trying again in {3 * attempt} s...")
            time.sleep(3 * attempt)


def main():
    # Some Windows consoles can't show every character in a title; print a ? instead of crashing
    sys.stdout.reconfigure(errors="replace")
    search = get_json("esearch.fcgi", {
        "db": "pubmed", "term": QUERY, "retmode": "json", "retmax": 10000,
        "datetype": "edat", "mindate": PARENT_SEARCH_END, "maxdate": SEARCH_END})
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
        writer = csv.DictWriter(f, fieldnames=[PLAIN[k] for k in rows[0]])
        writer.writeheader()
        writer.writerows({PLAIN[k]: v for k, v in r.items()} for r in rows)

    after = sum(r["after_parent_alerts"] for r in rows)
    print(f"{len(rows)} records added to PubMed from {PARENT_SEARCH_END} to {SEARCH_END}; "
          f"{after} after the parent's search alerts ended (~{ALERTS_END}).")
    print(f"Screened so far: {len(SCREENED)}. Wrote {OUT}")

    # Compare with the saved search used in the paper
    if SAVED.exists():
        with SAVED.open(encoding="utf-8") as f:
            saved = {r["PMID"] for r in csv.DictReader(f)}
        added = [r for r in rows if r["pmid"] not in saved]
        missing = saved - {r["pmid"] for r in rows}
        print(f"\nSaved search used in the paper ({SAVED}): {len(saved)} records.")
        if not added and not missing:
            print("This rerun found exactly the same records.")
        else:
            print(f"This rerun has {len(added)} record(s) the saved search did not, and is missing "
                  f"{len(missing)}. PubMed adds subject terms to records weeks after they first "
                  "appear, so older records can start matching the search later:")
            for r in added:
                print(f"  + PMID {r['pmid']}, added to PubMed {r['entrez_date']}: {r['title'][:70]}")
            for pmid in sorted(missing):
                print(f"  - PMID {pmid}")


if __name__ == "__main__":
    main()
