"""Write the sample sheet of one or more GEO series as a TSV: one row per GSM
with its title, the supplementary file names and every "characteristics"
field as its own column (donor, age, sex, tissue level and so on).

Usage: python workflow/scripts/geo_samples.py GSE201586 GSE168243 > data/metadata/geo_samples.tsv
"""
import re
import sys
import urllib.request

UA = {"User-Agent": "pnatlas-stage2"}


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def samples(gse):
    txt = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief")
    rows, cur = [], None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"series": gse, "gsm": line.split("=", 1)[1].strip(), "supp_files": []}
            rows.append(cur)
        elif cur is None:
            continue
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            val = line.split("=", 1)[1].strip()
            k, _, v = val.partition(":")
            key = re.sub(r"[^a-z0-9]+", "_", k.strip().lower()).strip("_")
            cur[key] = v.strip()
        elif line.startswith("!Sample_supplementary_file"):
            cur["supp_files"].append(line.split("=", 1)[1].strip().rsplit("/", 1)[-1])
    for r in rows:
        r["supp_files"] = ";".join(r["supp_files"])
    return rows


def main(accs):
    rows = [r for gse in accs for r in samples(gse)]
    cols = ["series", "gsm", "title"]
    cols += sorted({k for r in rows for k in r} - set(cols) - {"supp_files"}) + ["supp_files"]
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r.get(c, "")) for c in cols))


if __name__ == "__main__":
    main(sys.argv[1:])
