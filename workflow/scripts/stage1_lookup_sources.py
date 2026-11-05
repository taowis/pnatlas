"""Look up where a dataset's files and source paper are, for sources that the
GEO and non-GEO checks did not settle.

Subcommands:
  zenodo <record_id>        latest version, licence, file list with md5
  geo <GSE>                 summary, design, contributors, supplementary files
  epmc-accession <GSE>      Europe PMC articles mentioning the accession
  epmc-fulltext <PMCID> <regex>   lines of the open full text matching regex
  pubmed <PMID>             title, journal, date

Usage: python workflow/scripts/stage1_lookup_sources.py geo GSE189501
"""
import html
import json
import re
import sys
import urllib.parse
import urllib.request

UA = {"User-Agent": "pnatlas-stage1-check"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def zenodo(rec):
    d = json.loads(fetch(f"https://zenodo.org/api/records/{rec}"))
    m = d["metadata"]
    print(f"record\t{d['id']}\ttitle\t{m['title']}")
    print(f"doi\t{d.get('doi')}\tconcept_doi\t{d.get('conceptdoi', '')}")
    print(f"licence\t{m.get('license', {}).get('id')}\taccess\t{m.get('access_right')}"
          f"\tpublication_date\t{m.get('publication_date')}")
    for f in d.get("files", []):
        print(f"file\t{f['key']}\t{f['size']}\t{f.get('checksum')}")


def geo(acc):
    t = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief")
    keep = ("!Series_title", "!Series_summary", "!Series_overall_design", "!Series_contributor",
            "!Series_contact_institute", "!Series_supplementary_file", "!Series_relation",
            "!Series_pubmed_id", "!Series_sample_taxid", "!Series_submission_date")
    for line in t.splitlines():
        if line.startswith(keep):
            print(line)


def epmc_accession(acc):
    q = urllib.parse.urlencode({"query": acc, "format": "json", "pageSize": 10})
    d = json.loads(fetch(f"{EPMC}/search?{q}"))
    for h in d["resultList"]["result"]:
        print(f"{h.get('pmid', '')}\t{h.get('pubYear', '')}\t{h.get('journalTitle', '')}\t{h.get('title', '')}")


def epmc_fulltext(pmcid, pattern):
    t = fetch(f"{EPMC}/{pmcid}/fullTextXML")
    seen = set()
    for m in re.finditer(pattern, t):
        s = re.sub(r"<[^>]+>", " ", t[max(0, m.start() - 200):m.start() + 200])
        s = re.sub(r"\s+", " ", html.unescape(s)).strip()
        if s not in seen:
            seen.add(s)
            print("..", s)


def pubmed(pmid):
    d = json.loads(fetch(f"{EUTILS}/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"))
    s = d["result"][pmid]
    print(f"{pmid}\t{s.get('source')}\t{s.get('pubdate')}\t{s.get('title')}\t{s.get('elocationid', '')}")


def main(argv):
    cmd, args = argv[0], argv[1:]
    {"zenodo": zenodo, "geo": geo, "epmc-accession": epmc_accession,
     "epmc-fulltext": epmc_fulltext, "pubmed": pubmed}[cmd](*args)


if __name__ == "__main__":
    main(sys.argv[1:])
