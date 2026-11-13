# Static gene lookup site

`build_portal.py` reads the result tables of this repository and writes a
self-contained site: `index.html` (this directory's page, copied as is), the
JSON the page reads (a gene index, one shard of gene records per two-letter
symbol prefix, the state table, the site metadata) and copies of the result
tables offered for download. Nothing is computed in the browser and no server
is needed beyond a file host, so the site needs no maintenance once deployed.

```bash
python portal/build_portal.py            # writes portal/site/ (ignored by git)
cd portal/site && python -m http.server 8000   # then open http://localhost:8000
```

The pages workflow (`.github/workflows/pages.yml`) rebuilds and deploys the
site to GitHub Pages on every push to `main` that touches the portal, the
results or the metadata, once Pages is set to "GitHub Actions" in the
repository settings.

What the page shows per gene: specificity and log CPM in each of the 31
reference states, the compartment sums over the states mapped in tissue, the
correlation of the gene with each nerve state's fraction across the GTEx
tibial nerve samples, its colocalisation tests with the five neuropathy GWAS,
and its membership of the neuropathy gene classes and the PNPrank truth set.
The states page lists cells, donors, datasets, spatial mapping status, GTEx
mean fractions and enriched gene classes per state.
