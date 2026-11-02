# Results

| File | Content |
|---|---|
| `target_gene_statistics.csv` | Primary (Pnoc, Oprl1), secondary (Pcsk1/2, Cpe) and exploratory opioid tests for every contrast: n, mean normalized abundance, detection, log2FC, SE, 95% CI, p, family BH/Holm, transcriptome padj, result class |
| `target_gene_sample_values.csv` | Per-animal raw and normalized values for all target genes |
| `de/*_transcriptome.csv.gz` | Full transcriptome results for each contrast |
| `sensitivity_primary.csv` | NB GLM, leave-one-out, Mann–Whitney, sex-stratified and QC-exclusion checks |
| `exploratory_interactions.csv` | Diet × age (GSE225224, GSE320302) and sex × diet (GSE281091) terms |
| `qc/` | Per-sample QC (library size, size factors, PCA, similarity flags), target dispersions |
| `single_cell/` | GSE233163 pool × cell-type summaries and composition; GSE309905 context |
| `model_notes.json` | Annotation checks, duplicated IDs, Cook's distances, moderation priors |
| `feeding/` | Modules A–C: fasting/refeeding/HFD nodose models, panel statistics, sensitivities, concordance permutation, length-bias index |
| `environment_pip_freeze.txt` | Exact package versions used |

See `../REPORT.md` for interpretation.
