# T. pacificus Venom Transcriptomics — Pipeline Breakdown

---

## What This Pipeline Does

This pipeline analyzes RNA-seq data from 6 tissues in the Japanese flying squid (*Todarodes pacificus*) to identify and characterize venom genes, with a focus on comparing the two venom-producing glands — the **Posterior Salivary Gland (PSG)** and **Accessory Salivary Gland (ASG)**.

**The 6 tissues are:**
- PSG — Posterior Salivary Gland (venom)
- ASG — Accessory Salivary Gland (venom)
- BM — Buccal Mass
- SU_ms — Supraoesophageal Mass
- Arm
- Eso — Oesophagus

---

## Step-by-Step Breakdown

---

### Section 1 — Load Libraries
Load all R packages needed for the analysis. We use the `conflicted` package to prevent errors when multiple packages share the same function name (e.g. `dplyr::select` vs `AnnotationDbi::select`).

---

### Section 2 — Configuration
Define tissue names, file paths, and color schemes in one place so they stay consistent across every figure and table in the analysis.

---

### Section 3 — Load Kallisto Data
**Kallisto** was run separately (outside R) to quantify gene expression from raw RNA-seq reads. Here we import those results using `tximport`, which:
- Reads the `abundance.tsv` files from each tissue
- Converts transcript-level counts to gene-level counts using a transcript-to-gene mapping file (`tx2gene.tsv`)

---

### Section 4 — Normalize Counts
Raw read counts are not directly comparable between samples because sequencing depth differs. We use **DESeq2 size-factor normalization** to correct for this, then save the normalized counts for downstream use.

---

### Section 5 — VST Transformation & Fold-Change Filter
**Why VST?**
Raw counts have very high variance at low expression levels which distorts clustering and PCA. Variance Stabilizing Transformation (VST) compresses this variance so that all expression levels are treated more equally.

**Why the fold-change filter?**
Genes that are expressed at the same level in every tissue tell us nothing about tissue specialization. We keep only genes with at least a **2-fold difference** (1 unit on log2 scale) across any tissue — reducing ~185,000 genes to ~89,000 biologically informative ones.

---

### Section 6 — Load EggNOG Annotations
EggNOG-mapper was run separately to functionally annotate the Trinity transcriptome assembly. It assigns each gene:
- A functional description
- GO terms
- KEGG pathway IDs
- PFAM domain IDs

We load annotation files for all 6 tissues, clean the gene IDs to match our expression data, then deduplicate keeping the best annotation per gene.

**Result:** ~58,000 of 185,000 genes have functional annotations (~31%)

---

### Section 7 — Build GO & KEGG Databases
We restructure the EggNOG annotations into the format that `clusterProfiler` requires for enrichment analysis:
- A **term2gene** table — which genes belong to each GO/KEGG term
- A **term2name** table — human-readable names for each term

We split GO terms into the three ontologies:
- **BP** — Biological Process (what the gene does)
- **MF** — Molecular Function (how it does it)
- **CC** — Cellular Component (where in the cell)

---

### Section 8 — Ranked Gene Lists for GSEA
**Gene Set Enrichment Analysis (GSEA)** requires genes ranked by how strongly they are enriched in a tissue of interest.

For each tissue we calculate:
```
rank score = expression in target tissue − mean expression in all other tissues
```
A high positive score means a gene is specifically expressed in that tissue. A negative score means it is suppressed relative to other tissues.

---

### Section 9 — Run GSEA
GSEA asks: *"Do genes belonging to a particular biological pathway tend to cluster at the top or bottom of our ranked list?"*

If a pathway's genes are consistently at the top → that pathway is enriched in our target tissue.
If they cluster at the bottom → the pathway is depleted.

We run GSEA for:
- **PSG vs all other tissues**
- **ASG vs all other tissues**
- **PSG vs ASG directly**

Across GO BP, GO MF, GO CC, and KEGG Pathways.

---

### Section 10 — GSEA Visualizations
We generate three types of figures:
- **Dotplots** — show enrichment score and significance for top terms
- **NES Barplots** — show which pathways are enriched vs depleted (PSG = red, ASG = blue)
- **NES Heatmap** — shows enrichment patterns across all 6 tissues simultaneously

---

### Section 11 — Venom Gene Identification
We search the EggNOG annotations for genes likely to be venom-related using keyword lists:

| Category | Examples |
|---|---|
| **Venom toxins** | toxin, phospholipase, metalloprotease, lectin, kunitz, CRISP |
| **Venom processing** | furin, signal peptide, disulfide isomerase, SNARE, Rab, exocytosis |
| **Venom PFAM domains** | Known toxin structural domains (PF00087, PF00188 etc.) |

This gives us ~2,450 candidate venom-related genes.

---

### Section 12 — Venom Gene Expression Analysis
For each venom gene we calculate:
- **PSG enrichment** = PSG expression − mean of non-venom tissues
- **ASG enrichment** = ASG expression − mean of non-venom tissues
- **PSG vs ASG** = PSG − ASG directly

This lets us distinguish genes that are PSG-specific, ASG-specific, or shared between both glands.

---

### Section 13 — Venom Gene Clustering
We use **hierarchical clustering** (Ward's method) on the VST expression matrix of all venom genes to group them by expression pattern across tissues.

Genes are split into 5 clusters and each cluster is characterized by:
- Which tissue it is highest in
- What functional category dominates
- Mean PSG enrichment score

This reveals that venom genes fall into biologically distinct groups rather than being one uniform set.

---

### Section 14 — PCA of Venom Genes
**Principal Component Analysis (PCA)** reduces the 6-dimensional expression space (one dimension per tissue) down to 2 dimensions we can visualize.

- Each **point** is a venom gene, positioned by its expression pattern
- **Arrows** show where each tissue "pulls" — genes near the PSG arrow are PSG-enriched, genes near the ASG arrow are ASG-enriched
- **Point size** = how strongly the gene is enriched in PSG
- **Color** = which expression cluster the gene belongs to
- The **density inset** shows the distribution of all genes along PC1, confirming clear separation between PSG and ASG gene groups

PC1 (35% of variance) separates PSG from all other tissues. PC2 (24%) separates ASG from non-venom tissues.

---

### Section 15 — Final PCA Figure
Assembles the final publication-quality PCA plot with:
- Clean gene labels (all truncated and cryptic descriptions replaced with readable names)
- Deduplicated labels (no gene appears twice)
- ASCII-safe text (no Unicode arrows that cause PDF encoding warnings)
- Density inset in the bottom-left corner
- Saved as both PDF (for publication) and PNG (for presentations)

---

### Section 16 — Summary & Session Info
Prints a summary of all results and saves `session_info.txt` recording the exact R and package versions used — important for reproducibility.

---

## Output Files

| File | Description |
|---|---|
| `ALL_T_pacificus_Normalized_Genes.csv` | Normalized expression for all 185,320 genes |
| `eggnog_all_genes_annotated.csv` | Functional annotations for 58,407 genes |
| `eggnog_filtered_2fold_annotated.csv` | Annotations filtered to 2-fold variable genes |
| `custom_annotation_databases.RData` | GO and KEGG databases for GSEA |
| `GSEA_results.RData` | All GSEA result objects |
| `GSEA_BP/MF/CC/KEGG_PSG_vs_ASG.csv` | GSEA tables for each ontology |
| `venom_related_genes_annotated.csv` | All 2,450 candidate venom genes |
| `PSG_venom_genes_summary.csv` | PSG-enriched venom genes ranked |
| `venom_genes_clusters.csv` | Cluster assignments for all venom genes |
| `venom_genes_PCA_FINAL_clean.pdf/.png` | Final PCA figure |
| `GSEA_dotplot_*.pdf` | GSEA dotplots |
| `GSEA_NES_barplot_PSG_vs_ASG.pdf` | NES barplot |
| `GSEA_NES_heatmap_all_tissues.pdf` | Cross-tissue enrichment heatmap |
| `session_info.txt` | R session and package versions |
