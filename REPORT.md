# Hypercaloric feeding and *Pnoc*/*Oprl1* expression in peripheral ganglia — reanalysis report

Analysis date: 2 November 2026. All inputs are public GEO deposits retrieved on that date
(URLs, byte sizes and SHA-256 hashes in `metadata/download_manifest.json`).

## Summary

Across five testable diet contrasts in four independent mouse studies (DRG and nodose ganglia,
adult and early-life high-fat feeding), **neither *Pnoc* nor *Oprl1* showed a robust or reproducible
change in ganglion transcript abundance with hypercaloric feeding.**

- **One nominal signal:** DRG *Pnoc* was lower after 9 weeks of 45% HFD (age 12 wk; log2FC −3.16,
  95% CI −5.10 to −1.22; BH q = 0.014 in the prespecified Pnoc/Oprl1 family; mean normalized
  counts 39.2 vs 4.8). It is **fragile**: *Pnoc* is barely expressed in bulk DRG (mean ≈20
  normalized counts, zeros in 3/12 animals, dispersion 1.24), one control animal contributes most of
  the control mean, a rank test is not significant (Mann–Whitney p = 0.09), leave‑one‑out p reaches
  0.057, it is not transcriptome‑wide significant (padj = 0.20) and it does **not** reproduce after
  33 weeks of the same diet in the same study (log2FC −0.62, CI −2.57 to 1.32).
- ***Oprl1*** was stable everywhere it could be tested. All 95% CIs lie within roughly ±0.55 log2
  units (≈ ±45%); the largest point estimate is +0.28 (DRG, 12 wk; q = 0.21).
- **Non-testable:** in Nav1.8⁺ DRG TRAP (GSE198895) both genes have zero counts in all four
  libraries; in sympathetic scRNA-seq (GSE233163) *Pnoc* is absent from the deposited matrix and
  diet is completely confounded with sequencing run.
- **Localization context:** in sympathetic ganglia *Oprl1* is neuron-restricted (detected in
  0–41% of neurons per pool; ≤0.1% of every non-neuronal type). In HFD-fed vagus snRNA-seq
  (GSE309905, no control diet) *Oprl1* is detected in 8–14% of putative neuronal nuclei and *Pnoc*
  essentially not at all.

**Interpretation.** These public data do not support a consistent transcriptional up‑ or
down‑regulation of the N/OFQ precursor or the NOP receptor in somatosensory or vagal ganglia
during hypercaloric feeding. Bulk tissue data have limited sensitivity for *Pnoc* because of its low
abundance; a decrease confined to a small neuronal subpopulation, or post-transcriptional changes,
cannot be excluded. The DRG 12‑week *Pnoc* decrease is a hypothesis for targeted validation
(e.g. RNAscope or qPCR on independent animals), not an established effect.

## Data and eligibility

| Contrast | Series | Tissue / assay | Diet & timing | n C/H | Status |
|---|---|---|---|---|---|
| DRG_HFD_12wk | GSE225224 | lumbar DRG, bulk RNA-seq | 45% HFD vs 10% fat SD, 9 wk (age 12 wk) | 6/6 | inferential |
| DRG_HFD_36wk | GSE225224 | as above | 33 wk (age 36 wk) | 5/6 | inferential |
| DRG_Nav18_TRAP_WD | GSE198895 | Nav1.8⁺ DRG RiboTag TRAP | Western diet vs chow, 11 wk | 2/2 pools | exploratory; targets not detected |
| Nodose_HFD_4wk | GSE281091 | nodose, bulk RNA-seq | 60% HFD vs non-fasted chow, 4 wk; 3 M + 3 F per arm | 6/6 | inferential (normalized data only) |
| Nodose_earlyHFD_P21 | GSE320302 | nodose, bulk RNA-seq | reared on HFD, at weaning | 8/8 | inferential; litter not deposited |
| Nodose_earlyHFD_P35 | GSE320302 | as above | 2 wk after weaning to chow (persistence) | 8/7 | inferential; litter not deposited |
| Symp_obese_scRNA | GSE233163 | SCG + stellate, 10x scRNA-seq | 60% HFD vs chow, 12 wk | pools | descriptive only |
| — | GSE309905 | vagus/nodose snRNA-seq | HFD only (Pten WT/KO) | — | context only |

Exclusions: GSE225224 sciatic nerve (not a ganglion); GSE281091 fasted and fed-state chow arms
(not a diet contrast); GSE309905 (no control diet). The Dryad DRG single-cell lead
(doi:10.5061/dryad.9s4mw6mm5) still returns 404 from doi.org and the Dryad API. Full sample-level
metadata: `metadata/sample_manifest.csv`; frozen contrast list: `metadata/contrast_inventory.csv`
(committed before any target effect was estimated).

## Methods

- **Count data (GSE225224 featureCounts; GSE198895 htseq-count; GSE320302 RSEM expected counts,
  rounded):** PyDESeq2 0.5.2, design `~group` (age × diet cells), median-of-ratios size factors and
  dispersions estimated on the full transcriptome after removing all-zero features only; Wald tests
  of the prespecified contrasts; unshrunk MLE log2 fold changes with nominal 95% Wald CIs; Cook's
  outlier filtering and independent filtering at default for transcriptome-wide padj.
- **GSE281091:** only DESeq2 median-ratio normalized values (Partek, pre-filtered) are deposited and
  they cannot be converted back to integers, so a limma-trend–style empirical-Bayes moderated linear
  model on log2(normalized + 1) with `~ sex + diet` was used (prior df 3.2; residual df 9).
  Pseudocount 0.5 sensitivity is in `results/target_gene_statistics.csv`.
- **Multiplicity:** Benjamini–Hochberg across the prespecified primary family (Pnoc, Oprl1 × the
  contrasts with an estimate = 10 tests); Holm as conservative sensitivity. Pcsk1/Pcsk2/Cpe and the
  opioid panel (Penk, Pdyn, Pomc, Oprm1, Oprd1, Oprk1) are separate secondary/exploratory families.
- **Single cell (GSE233163):** dense UMI matrix streamed to pseudobulk sums by library × deposited
  cell type; libraries merged to site × run pools (ganglia were pooled per condition; 12k/14k
  loadings are not independent animals). Coverage rule (set before inspecting targets): ≥20 cells
  per pool × cell type. No diet test: lean pools come from runs N61/N66 and obese pools only from N68.
- **Robustness (defined after the primary results were seen; reported, not used for selection):**
  single-gene NB GLM with fixed MAP dispersion (reproduces PyDESeq2 exactly), leave-one-animal-out
  refits, Mann–Whitney U, diet × age interactions, and for GSE281091 sex-stratified fits and a fit
  excluding three QC-flagged samples (`results/sensitivity_primary.csv`, `results/exploratory_interactions.csv`).

## Primary results (`results/target_gene_statistics.csv`)

| Contrast | *Pnoc* log2FC (95% CI) | q | *Oprl1* log2FC (95% CI) | q |
|---|---|---|---|---|
| DRG HFD, 12 wk | −3.16 (−5.10, −1.22) | 0.014 | +0.28 (0.01, 0.54) | 0.21 |
| DRG HFD, 36 wk | −0.62 (−2.57, 1.32) | 0.59 | −0.14 (−0.41, 0.14) | 0.51 |
| Nav1.8⁺ TRAP, WD | not detected | — | not detected | — |
| Nodose HFD, 4 wk | +0.32 (−0.49, 1.13) | 0.51 | −0.08 (−0.21, 0.06) | 0.51 |
| Nodose early-life HFD, P21 | −0.39 (−1.22, 0.44)ᴸ | 0.51 | +0.08 (−0.24, 0.39) | 0.64 |
| Nodose early-life HFD, P35 | +0.48 (−0.53, 1.49)ᴸ | 0.51 | −0.20 (−0.53, 0.13) | 0.51 |

ᴸ low abundance (mean normalized count < 10). Positive values = higher under hypercaloric diet.
Non-significance is not evidence of no change; CI widths show what each study could exclude.

Secondary genes: no *Pcsk1*, *Pcsk2* or *Cpe* contrast reached q < 0.4. Exploratory opioid panel:
nothing reached q < 0.3 (lowest: nodose *Oprm1* −1.69, *Pomc* +0.54, both nominal p ≈ 0.02, q = 0.31).

Transcriptome-wide context (genes with padj < 0.05): DRG 12 wk 23, DRG 36 wk 0, TRAP 189,
nodose adult 0, nodose P21 33, nodose P35 0 — i.e. most diet contrasts produce little bulk
transcriptional signal in these tissues.

## Quality control and caveats

- **GSE281091:** PC1 (62% of top-variance-gene variation, dominated by Gm/non-coding features)
  separates the three HFD males from all other samples, so the male diet effect is confounded with a
  technical axis. Excluding them or restricting to females does not change the *Pnoc*/*Oprl1*
  conclusion (female-only *Pnoc* −0.47, p = 0.23; *Oprl1* −0.04, p = 0.31). The nominal sex × diet
  interaction for *Pnoc* (p = 0.044, unmoderated) is therefore not interpretable.
- **GSE225224:** one 12-wk HFD sample (GSM7042132) has lower similarity and fewer detected genes;
  it was retained (removal would not be based on a pre-stated rule). No spinal-cord contamination
  signature (*Gad1/2*, *Slc32a1*, *Pax2*, *Mog*) tracks DRG *Pnoc* across samples. Sequencing batch
  is not deposited. Diet × age interaction for *Pnoc*: +2.53, p = 0.07 (exploratory).
- **GSE320302:** litter membership is not deposited; pups from the same dam may not be independent,
  so CIs may be too narrow. Four samples flagged by QC were retained.
- **GSE198895:** n = 2 pools per arm; *Pnoc*/*Oprl1* features exist in the count files but are 0
  in every library.
- **GSE233163:** *Pnoc*, *Oprd1*, *Oprk1* and *Pdyn* are absent from the deposited (pre-filtered)
  17,952-gene matrix; only 406 neurons overall.
- Process note: while checking whether GSE281091 normalized values could be converted back to raw
  counts, the per-sample *Pnoc*/*Oprl1* values were printed before the contrast inventory was frozen.
  No effect estimate was computed at that point, and the inventory follows the pre-existing plan.
- Analyses ran in Python 3.13 with the pinned packages in `requirements.txt` (full freeze:
  `results/environment_pip_freeze.txt`); the plan specified 3.12.

## Figures (legends in `figures/LEGENDS.md`; source data in `figures/source_data/`)

- **Fig. 1** *Pnoc* and *Oprl1* diet effects by ganglion, exposure and assay (forest plots).
- **Fig. 2** Biological-replicate normalized abundance for the five testable contrasts.
- **Fig. 3** Single-cell localization context (sympathetic and vagal).
- **Fig. S1** Sample QC (PCA, library sizes). **Fig. S2** Pcsk1/Pcsk2/Cpe. **Fig. S3** Opioid panel.

Fonts: Nimbus Sans converted from the URW OTF (CFF) files to TrueType outlines at build time;
`pdffonts` shows only embedded CID TrueType Nimbus Sans subsets, with no Poppler mismatch warning.

## Reproduce

```bash
bash scripts/run_all.sh
```
