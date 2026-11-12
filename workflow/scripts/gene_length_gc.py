"""Per-gene length and GC content from an Ensembl GTF and the matching cDNA
FASTA, for the length and GC matched background of the stage 5
enrichment (plan section 6). Length is the number of genomic bases
covered by the union of the gene's exons; GC is the fraction of G and C
over the concatenated sequences of the gene's transcripts in the cDNA
file (every transcript once, so shared exons weigh by the number of
transcripts that carry them). Only chromosomes 1 to 22, X and Y are kept.

Writes a gzipped table: ensembl_gene_id, gene_name, gene_biotype,
chromosome, exon_length, gc_percent, transcripts.
Usage: python workflow/scripts/gene_length_gc.py --gtf <Homo_sapiens.GRCh38.105.gtf.gz> \
           --cdna <Homo_sapiens.GRCh38.cdna.all.fa.gz> --out data/metadata/ensembl105_gene_length_gc.tsv.gz
"""
import argparse
import gzip
import re
from collections import defaultdict

import pandas as pd

CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def merged_length(intervals):
    intervals.sort()
    total, cur_s, cur_e = 0, None, None
    for s, e in intervals:
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--cdna", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    attr_re = re.compile(r'(\S+) "([^"]*)"')
    exons = defaultdict(list)
    info = {}
    with gzip.open(args.gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] not in CHROMS or p[2] not in ("gene", "exon"):
                continue
            a = dict(attr_re.findall(p[8]))
            gid = a["gene_id"]
            if p[2] == "gene":
                info[gid] = (a.get("gene_name", ""), a.get("gene_biotype", ""), p[0])
            else:
                exons[gid].append((int(p[3]), int(p[4])))
    gc = defaultdict(lambda: [0, 0, 0])   # gc bases, all bases, transcripts
    with gzip.open(args.cdna, "rt") as f:
        gid, seq = None, []

        def flush():
            if gid is not None and gid in info:
                s = "".join(seq).upper()
                gc[gid][0] += s.count("G") + s.count("C")
                gc[gid][1] += len(s)
                gc[gid][2] += 1

        for line in f:
            if line.startswith(">"):
                flush()
                m = re.search(r"gene:(ENSG\d+)", line)
                gid, seq = (m.group(1) if m else None), []
            else:
                seq.append(line.strip())
        flush()
    rows = []
    for g, (name, biotype, chrom) in info.items():
        b = gc.get(g)
        rows.append([g, name, biotype, chrom, merged_length(exons[g]) if g in exons else 0,
                     round(100 * b[0] / b[1], 2) if b and b[1] else float("nan"), b[2] if b else 0])
    df = pd.DataFrame(rows, columns=["ensembl_gene_id", "gene_name", "gene_biotype", "chromosome", "exon_length",
                                     "gc_percent", "transcripts"]).sort_values("ensembl_gene_id")
    df.to_csv(args.out, sep="\t", index=False)
    pc = df[df["gene_biotype"] == "protein_coding"]
    print(f"{len(df)} genes ({len(pc)} protein coding, {int(pc['gc_percent'].notna().sum())} with cDNA GC) written to {args.out}",
          flush=True)


if __name__ == "__main__":
    main()
