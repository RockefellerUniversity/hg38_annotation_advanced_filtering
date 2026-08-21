# hg38 annotation: advanced filtering for nuclear RNA-seq

Confident gene sets for nuclear RNA-seq that is counted over whole gene bodies. Apply one of these
sets to your count matrix **before** the DESeq2 step, so differential expression runs only on genes
we trust.

## Where the filters are

* `output/genes_*.txt` — one gene set per file, a plain list of Ensembl gene IDs (one per line).
* `output/gene_confidence_table.csv` — the full per-gene table; each filter is a logical column
  (`keep`, `keep_high`, ...), alongside all the evidence it is based on.
* `reports/GeneSet_Filtering.html` — how the sets are built and what each column means.
* `reports/AnnotationComparison.html` — BRC current annotation vs current GENCODE.

The filters, least to most strict (each a subset of the one before):

| Filter | Genes | What it is |
|---|---|---|
| `keep` | ~32.6k | current + wanted biotype (protein-coding, lncRNA) |
| `keep_high_medium` | ~32.2k | `keep` with some support evidence |
| `keep_high` | ~24.2k | `keep` at the high-confidence tier |
| `*_no_nested` | — | any of the above, with genes nested inside another kept gene removed |

## Applying a filter before DESeq2

Subset the count matrix to the chosen gene set, then build the DESeq dataset:

```r
genes  <- readLines("output/genes_keep_high.txt")
counts <- counts[rownames(counts) %in% genes, ]          # rownames are Ensembl gene IDs
dds    <- DESeq2::DESeqDataSetFromMatrix(counts, colData, design = ~ condition)
dds    <- DESeq2::DESeq(dds)
```

## Recommended filter and caveat

For a protein-coding differential-expression analysis use **`keep_high`**
(`output/genes_keep_high.txt`).

Caveat: `keep_high` is evidence-only, so it can still include a gene whose body sits inside another
kept gene. If that overlap ambiguity matters, use `keep_high_no_nested` instead — but note it drops
~12% of `keep_high`, and some of those are legitimate nested or bidirectional gene pairs. `keep_high`
is also protein-coding-focused; if you need lncRNAs, start from `keep_high_medium` and treat the
lncRNAs case by case (their evidence is thinner and they drive most of the overlap).
