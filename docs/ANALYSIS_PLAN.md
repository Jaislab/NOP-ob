# Hypercaloric feeding and the peripheral-ganglion N/OFQ–NOP axis

## Status

The uploaded workbook is a catalogue, not an expression dataset. Dataset contents,
replication, diet assignments and publications have not yet been independently
verified. Initial official GEO requests were denied by network policy. No gene
expression estimates, differential-expression tests or biological data figures
have been produced. Files named STYLE_PROOF_NOT_DATA are typography checks only.

## Primary question and estimands

Estimate changes in **Pnoc** and **Oprl1** transcript abundance associated with
hypercaloric feeding relative to each study's appropriate control diet. Report
positive log2 fold changes as higher abundance under hypercaloric feeding.
Distinguish ongoing adult feeding, lactational exposure and persistence after
return to chow. The analysis does not measure N/OFQ peptide abundance, release,
receptor protein, receptor activity or a causal transcriptional mechanism.

## Candidate contrasts (from the supplied catalogue; require verification)

| GEO series | Comparison | Important design checks |
|---|---|---|
| GSE225224 | DRG HFD vs SD at 12 and 36 weeks | Exclude sciatic nerve; verify independent animals, batch, ages and duration. If identifiable, fit diet × age and extract each age contrast. |
| GSE198895 | Nav1.8+ DRG Western diet vs chow | TRAP measures ribosome-associated RNA. Establish pooling and biological units; reported n=2/arm makes inference exploratory. |
| GSE281091 | Nodose HFD vs non-fasted chow | Exclude fasted/refed arms from the primary contrast; adjust for sex. A sex interaction is exploratory and requires identifiable replication. |
| GSE320302 | Nodose early-life HFD vs control at P21 and P35 | Establish dam/litter membership and independence. P35 is post-exposure persistence, not concurrent feeding. Use litter as the exposure unit or account for it with a supported model. |
| GSE233163 | Sympathetic HFD vs control within cell type and ganglion site | Establish pool/library/animal identity and batch overlap. Do not count cells, technical libraries or ganglia from the same animal as independent animals. If diet is confounded with experiment, do not claim an adjusted diet effect. |

GSE309905 lacks a reported control-diet arm and cannot estimate a diet effect.
Other candidate DRG datasets require an availability check before inclusion.
Genotype-only, STZ-only, non-ganglion and in-vitro lipid comparisons are outside
the primary question. No between-study chow-versus-HFD comparison will substitute
for a within-study diet contrast.

## Acquisition and eligibility

Retrieve official series/sample metadata and supplementary file listings; retain
URLs, download times and SHA-256 hashes. Review sample units and all relevant
characteristics before choosing analysis inputs. Prefer processed raw gene-count
matrices. Audit gene identifiers, reference annotation, strandedness when stated,
count semantics, duplicated identifiers, missing targets and sample-name mapping.
Do not round TPM/FPKM to fit a count model. If only normalized values are available,
use a method appropriate to those values or flag that count-based reanalysis needs
raw-read processing. Assess read volume and resources before such processing.

## Bulk/TRAP inference

Fit the negative-binomial model to the full usable transcriptome, not to the two
target genes alone; estimate library size factors and dispersions with PyDESeq2
0.5.2 where its model supports the verified design. Remove only wholly unexpressed
features initially; document any additional filtering. Inspect library sizes,
sample similarity/PCA, metadata consistency and outlier diagnostics. Do not remove
an outlier solely to improve a target-gene P value. Retain zeros and mark targets
that cannot be tested or have unstable estimates.

Report biological n per arm, mean normalized abundance, proportion of samples with
nonzero raw counts, log2 fold change, standard error, nominal 95% confidence interval,
raw P value and adjusted P value. Distinguish missing, zero-count, non-testable and
non-significant results. Non-significance alone is not evidence of no regulation.

For the primary family, apply Benjamini–Hochberg correction jointly to the Pnoc
and Oprl1 tests across all eligible, prespecified diet/time/cell-type contrasts;
freeze that inventory from metadata and annotation before looking at target
effects. Account for dependence in interpretation and provide a conservative
sensitivity correction if the final contrast structure warrants it. Confidence
intervals are nominal, not simultaneous. Also retain transcriptome-wide adjusted
P values as context. Clearly label exploratory interactions and secondary tests.

## Single-cell inference and composition

Inspect deposited cell annotations and quality-control information; preserve
sample identities. Plot expression/detection descriptively by cell type and
independent biological replicate. For differential expression, aggregate raw
counts by biological unit × cell type × ganglion site, using a prespecified cell
coverage rule chosen without consulting Pnoc/Oprl1 effects. A minimum of independent
replication in both diet arms and an identifiable design is required. Technical
libraries from a shared pool must be combined. If donor/pool identities cannot be
recovered, restrict conclusions to descriptive expression. No cell-level diet
P values. Assess cell-type proportions separately; bulk tissue cannot establish
within-cell transcriptional regulation from composition markers alone.

## Secondary genes and synthesis

Peptide-processing genes **Pcsk1**, **Pcsk2** and **Cpe** provide secondary context;
they are not specific readouts of NOP activity. Related opioid ligands/receptors
can be added as a separately declared exploratory panel. Avoid an invented
"NOP activation score" from nonspecific signaling genes. Keep bulk RNA and TRAP
estimands separate, and do not pool effects across different ganglia or developmental
exposures simply to obtain a single significant result.

## Figures and source data

- Main effect summary: Pnoc and Oprl1 forest plots, grouped by ganglion, exposure
  and assay, showing log2 fold changes with nominal 95% confidence intervals and
  sample sizes. No pooled effect unless scientifically and statistically justified.
- Sample-level panels: normalized abundance with individual biological-unit
  points. Control = blue **#0072B2** (circle); hypercaloric = pink **#CC79A7**
  (triangle), from the supplied PDF's Okabe–Ito palette. Use a stated transform,
  if required; do not use arbitrary pseudocounts to infer fold changes.
- Single-cell panels if eligible: replicate-aware cell-type expression/detection
  and pseudobulk diet effects, with cell coverage disclosed.
- Supporting panels: sample QC and secondary-gene results where informative.

Use Nimbus Sans, 7-pt text at final size (6.5-pt ticks/legends), lowercase bold
panel letters, legible axes and 0.6-pt axis rules. Use 89-mm single-column or
183-mm double-column canvases as appropriate. Export vector PDF with embedded
fonts, editable SVG and 600-dpi PNG. Preserve physical sizes, expose individual
replicates, avoid redundant significance stars, and provide exact statistics
and plotting values in source-data tables. This is a Nature-oriented production
specification, not a claim of journal acceptance or compliance with every title's
submission requirements.

## Execution handoff

See [HANDOFF.md](../HANDOFF.md) for current status, portable setup commands,
known limitations and the remaining execution checklist. All candidate study
metadata in this document come from the supplied catalogue and remain unverified.
