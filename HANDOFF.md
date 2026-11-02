# Execution handoff

## User requirements

Execute a reanalysis of publicly available peripheral-ganglion transcriptomes,
focusing on changes in **Pnoc** and **Oprl1** during hypercaloric feeding relative
to matched control diets. Examine other N/OFQ–NOP-related genes only as secondary
context. Deliver publication-quality figures with Nimbus Sans and the specified
blue/pink palette, reproducible code, source-data tables and a biological report.

The user approved the diet-focused plan and explicitly asked for execution.
They subsequently requested that the plans be pushed here so another AI can
execute them. Do not stop at another plan or a package installation.

## Work actually completed

- Inspected the supplied workbook and extracted all five sheets into CSV.
- Wrote the candidate dataset audit and detailed analysis specification.
- Installed PyDESeq2 0.5.2, requests 2.32.5, anndata 0.12.6 and statsmodels 0.14.6
  in an isolated workspace environment; imported those libraries successfully.
- Confirmed Nimbus Sans is available in the original cloud image and prepared a
  Matplotlib style helper. The generated specimen contains no biological data.
- Attempted official GEO metadata and supplementary-file requests, including a
  retry after the environment configuration changed to unrestricted networking.
- Verified GitHub access and prepared this repository handoff.

## Work not completed

**No official dataset/sample records or count matrices were acquired. No sample
manifest was verified. No expression model, hypothesis test, effect estimate,
biological conclusion or data figure exists.** There is no completed
differential-expression pipeline. The included Python scripts are acquisition,
catalogue-extraction and styling helpers only.

## Next execution steps

1. On a machine with working NCBI access, install `requirements.txt`, then run
   `scripts/acquire_geo.py`. If a series does not resolve, record that finding;
   do not assume the supplied accession or publication metadata is correct.
2. Inspect the official series/sample records and linked papers. Reconcile diet,
   tissue, age, sex, litter, batch, pooling and technical versus biological
   replication. Record all exclusions and URLs. Freeze an eligible contrast
   inventory before inspecting Pnoc/Oprl1 effect estimates.
3. Download deposited **raw gene counts** where available, retain SHA-256 hashes,
   and identify the reference annotation. The helper does not download count
   matrices automatically. Never treat the workbook as expression data.
4. Build a sample manifest with at least: series, sample accession, matrix column,
   tissue/site, assay, diet, exposure duration, collection age, sex, batch,
   biological unit/pool, technical replicate, dam/litter and inclusion reason.
   Unknown fields must remain explicitly unknown rather than guessed.
5. Run QC and implement the supported models described in
   [the analysis plan](docs/ANALYSIS_PLAN.md). Estimate normalization and dispersion
   from the usable transcriptome. Handle raw counts and already normalized data
   with appropriate, different methods. Resolve ambiguous gene mappings.
6. Analyze Pnoc/Oprl1 diet effects with biological replication. Keep adult HFD,
   Western diet/TRAP and lactational exposure separate. Do not infer a diet effect
   from a genotype-only comparison or compare unrelated control and HFD studies.
7. For single-cell data, recover donor/pool identity and experiment/batch
   structure before inference. Aggregate technical libraries from the same pool;
   cells are not independent animals. A confounded or unreplicated design supports
   descriptive results only. Document the coverage rule without selecting cells
   or contrasts for a desired target-gene effect.
8. Produce the figures and complete statistics/source-data tables, then inspect
   every panel at final physical size. Check font embedding and resolve the known
   Matplotlib CFF-font warning described below.
9. Write a concise report separating supported changes, uncertain or low-count
   estimates, non-testable targets and inaccessible studies. Save scripts,
   versions, provenance and commands required to reproduce each output.

## Dataset priority

| Priority | Catalogue candidate | Intended comparison |
|---|---|---|
| 1 | GSE225224 | DRG HFD vs SD at each reported age; omit sciatic nerve |
| 1 | GSE198895 | Nav1.8+ DRG Western-diet vs chow TRAP; reported n=2/arm, exploratory |
| 1 | GSE281091 | Nodose HFD vs non-fasted chow, sex-adjusted if design supports it |
| 1 | GSE320302 | Nodose early-life HFD vs control at P21 and P35; establish litter units |
| 2 | GSE233163 | Sympathetic diet effects within cell type/site, subject to replication and batch checks |
| Context | GSE309905 | HFD-only vagal expression localization; no reported diet contrast |
| Availability lead | Dryad DOI 10.5061/dryad.9s4mw6mm5 | Recheck the candidate DRG single-cell diet dataset |

All rows above are catalogue leads, not independently validated study designs.
Other inaccessible leads are preserved in `metadata/catalogue_datasets.csv`.

## Practical limitations of the helpers

- `acquire_geo.py` was prepared but has not completed a successful GEO retrieval.
  It validates the series marker in SOFT files, hashes successful downloads and
  returns a nonzero status if a request fails. It downloads metadata and listings,
  not expression files, and does not prove sample independence.
- `prepare_catalogue.py` can be rerun against the included workbook. It intentionally
  retains the supplied claims as unverified catalogue content.
- `figure_style.py` makes a clearly labelled **STYLE_PROOF_NOT_DATA** specimen.
  It is not an experimental result. Nimbus Sans must be installed; the script
  raises an error rather than silently substituting another family.
- Original font tests embedded Nimbus Sans, but Poppler emitted
  `Mismatch between font type and embedded font file` for Matplotlib's CFF/OTF
  export. Resolve that export issue before publication (for example, use a
  properly supported Nimbus Sans TrueType font or a verified vector renderer),
  then recheck with `pdffonts`, `pdftotext` and visual inspection. Do not claim
  the present helper alone guarantees publication-ready PDFs.
- The pinned direct dependency list is not a complete dependency lock, and no
  clean-machine installation or end-to-end reanalysis has been validated.

## Access blocker in the original environment

See [ACCESS_STATUS.md](docs/ACCESS_STATUS.md). Use the inherited proxy and valid TLS
verification where required. A new executor may have normal NCBI connectivity;
the failed original environment is not evidence that these datasets are private
or unavailable. Do not spend the study inventing workarounds for a broken cloud
instance when an authorized, functioning analysis machine is available.
