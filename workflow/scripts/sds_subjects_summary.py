"""Print donor count, sex and age range from a SPARC SDS subjects.xlsx.

Reads the xlsx with the standard library only (no pandas/openpyxl needed).
Usage: python workflow/scripts/sds_subjects_summary.py subjects.xlsx
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group()
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_rows(path):
    with zipfile.ZipFile(path) as z:
        strings = ["".join(t.text or "" for t in si.iter(NS + "t"))
                   for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for r in sheet.iter(NS + "row"):
        row = {}
        for c in r.iter(NS + "c"):
            v = c.find(NS + "v")
            if v is None:
                continue
            val = strings[int(v.text)] if c.get("t") == "s" else v.text
            row[col_index(c.get("r"))] = val
        if row:
            rows.append([row.get(i, "") for i in range(max(row) + 1)])
    return rows


def main(path):
    rows = read_rows(path)
    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}
    subs = [r for r in rows[1:] if r and r[0].startswith("sub-")]
    sexes = [r[idx["sex"]] for r in subs]
    ages = [int(re.match(r"\d+", r[idx["age"]]).group()) for r in subs]
    print(f"donors\t{len(subs)}")
    print(f"male\t{sexes.count('Male')}")
    print(f"female\t{sexes.count('Female')}")
    print(f"age_range\t{min(ages)} to {max(ages)} years")


if __name__ == "__main__":
    main(sys.argv[1])
