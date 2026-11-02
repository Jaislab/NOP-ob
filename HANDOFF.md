# Execution status

## Completed (2 November 2026)

- Official GEO records and supplementary files retrieved for all six catalogue series; SHA-256
  hashes in `metadata/download_manifest.json`. Original access blocker resolved (docs/ACCESS_STATUS.md).
- Sample manifest built from GEO SOFT records; contrast inventory frozen and committed before any
  target effect estimate (`git log -- metadata/contrast_inventory.csv`).
- Transcriptome-wide models: PyDESeq2 0.5.2 for GSE225224, GSE198895, GSE320302; moderated linear
  model for GSE281091 (raw counts not deposited, not reconstructable).
- Primary/secondary/exploratory target statistics with family-wise BH and Holm; robustness checks;
  QC; replicate-aware single-cell descriptions (GSE233163, GSE309905).
- Six figures (3 main, 3 supplementary) with legends and source data; inspected at final size;
  fonts verified with `pdffonts`/`pdftotext` (TrueType Nimbus Sans only; CFF mismatch resolved by
  converting the URW OTF outlines to TrueType at build time).
- Biological report: `REPORT.md`.

## Open limitations / possible next steps

- GSE281091 raw counts: reprocessing FASTQs from SRA would allow a count-based model and could
  clarify the technical PC1 axis affecting the three HFD-male samples.
- GSE320302 litter membership is not deposited; ask the authors to enable litter-aware models.
- GSE233163: diet confounded with run and Pnoc absent from the deposited matrix; the 2.7 GB Seurat
  object might contain the unfiltered gene set (not downloaded).
- Dryad DRG scRNA-seq lead (doi:10.5061/dryad.9s4mw6mm5) still 404; contact the Menichella lab.
- The DRG 12-wk Pnoc decrease should be validated in independent animals before being reported
  as an effect.
- Clean-machine installation validated only in this container (Python 3.13; pinned direct packages).
