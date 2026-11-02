# Feeding regulation in the vagal (nodose) ganglion — modules A–C

Analysis date: 2 November 2026. Gene panels, contrasts and multiplicity families were committed
before estimation (`metadata/feeding_gene_panels.csv`, `metadata/feeding_contrast_inventory.csv`,
`metadata/feeding_analysis_rules.md`; commit 0e2c7c9). Results: `results/feeding/`.

## Summary

1. **Fasting reshapes the nodose transcriptome; 4 weeks of HFD does not mimic it.** An 18-h fast
   changed 395 genes (padj < 0.05; 1,212–1,397 once a technical artifact is handled). The HFD
   signature does not resemble fasting: Spearman ρ = 0.09 across 19,711 genes, against a
   shared-control permutation null with mean 0.27 (P = 0.59). Within these data, the hypothesis that
   diet-induced obesity locks vagal afferents into a "fasted" state is not supported.
2. ***Pnoc* and *Oprl1* are not clearly state-sensitive.** With fasting, *Pnoc* −0.52 log2
   (95% CI −1.25 to 0.22) and *Oprl1* +0.11 (−0.02 to 0.23), q = 0.21 and 0.18. Females alone show
   the same directions more strongly (*Pnoc* −1.01, P = 0.043; *Oprl1* +0.19, P = 0.047; q = 0.09).
   This is a weak hint at reduced N/OFQ precursor and higher receptor during energy deficit, not a
   finding.
3. **Robust adult-HFD changes in nodose feeding genes:**
   - ***Uts2b*:** +1.07 log2 (CI 0.75–1.39), q = 3 × 10⁻⁵
   - ***Vip*:** +0.56 (0.31–0.81), q = 0.004
   - ***Lepr*:** −0.66 (−0.96 to −0.35), q = 0.004

   All three hold in both sensitivity analyses: female-only, and adjusted for gene-length bias
   (length-adjusted q = 0.0001, 0.018 and 0.019; *Lepr* is weaker in females alone, −0.24,
   P = 0.048). *Vip* and *Uts2b* are markers that published nodose single-cell atlases associate with
   a gut-mucosal vagal afferent subtype; lower *Lepr* fits reduced vagal leptin sensing in obesity.
   None of the three changed after early-life HFD (P21 or P35).
4. **Fasting-panel changes in feeding-relevant receptors** became clear once the artifact was
   handled: *Htr3a* ↑ (+0.32; length-adjusted q = 0.012), *Oxtr* ↓ (−0.37; female-only q = 0.036)
   and *Cnr1* ↑ (+0.21; length-adjusted q = 0.08). The CB1 increase matches the classic rat-nodose
   fasting result; *Cartpt*, *Mchr1* and *Npy2r* showed no fasting switch in mouse bulk ganglia.
5. **20-min refeeding** after an 18-h fast produced no detectable transcriptional change
   (0 genes at padj < 0.05; immediate-early genes all q ≥ 0.25). This is expected: immediate-early
   mRNAs peak 30–60 min after a stimulus and are diluted in bulk ganglion RNA.
6. **Early-life HFD (GSE320302):** no satiety-panel gene reached q < 0.05 at P21 or P35. The lowest
   were *Ffar3* −0.84 at P35 (P = 0.015, q = 0.08) and *Npy1r* −0.71 at P35 (P = 0.046).

## Important technical finding: a gene-length artifact in GSE281091

PC1 of the 24-sample nodose dataset (47% of top-gene variance) is a gene-length axis. Genes that
fall along it have a median length of 90 kb, versus 7 kb for genes that rise — the typical
signature of degraded RNA or 3′-coverage bias. A per-sample length-bias index puts the three
HFD-male samples at the extreme (−0.34 to −0.38), with one fed female and one fasted male partly
affected (Fig. S4).

Because the artifact coincides with HFD in males, several apparent HFD effects in the prespecified
model are artifacts. *Cartpt* ↑, *Cnr1* ↓, *P2rx3* ↓, *Insr* ↓ and *Piezo1* ↓ all correlate
|ρ| ≥ 0.8 with PC1 and disappear in females alone. Under either sensitivity, the HFD contrast has
0 transcriptome-wide significant genes. Only HFD effects that survive both sensitivities are
reported as robust above.

The two sensitivity models were added after this QC finding, not prespecified; they are reported
alongside the prespecified results (`results/feeding/feeding_panel_sensitivity*.csv`). The same
artifact does not alter the *Pnoc*/*Oprl1* conclusions of the primary report.

## Exploratory: what the fasting signature contains (post hoc reading, not a formal enrichment test)

- **Up:**
  - Glucocorticoid-responsive genes: *Fkbp5* (+0.61, padj 0.008), *Tsc22d3*/GILZ, *Hif3a* (+1.88),
    *Zbtb16*, *Klf9*, *Klf15*, *Mt1*/*Mt2*
  - *Plin4* (+3.7, a lipid-droplet gene)
  - Cold-inducible RNA-binding proteins *Cirbp*/*Rbm3*, consistent with the fall in body temperature
    mice show during fasting
- **Down:**
  - Injury/activity-associated genes: *Sprr1a* (−4.8), *Atf3*, *Ecel1*, *Gpr151*, *Arc*
  - Neuropeptides and receptors: *Ucn*, *Nts*, *Grp*, *Cckbr*
  - MHC-II genes: *H2-Aa*, *H2-Ab1*, *Cd74*
- **Interpretation:** an energy-deficit/corticosterone-driven program plus reduced
  activity-dependent gene expression in vagal afferents when the gut is empty.

## Methods (brief)

- **Combining the two GSE281091 files.** The deposited normalized files (3A: fasted, refed, fed;
  3B: HFD, fed) share the six fed samples. 3A was rescaled onto 3B with a single factor (0.989;
  maximum deviation 0.008 log2 on the shared samples), giving 21,621 common genes.
- **Model.** A moderated linear model on log2(normalized + 1), `~ group + sex`, 24 samples
  (residual df 19, prior df 3.6). It is the same method as the primary nodose analysis, with
  variance estimated across all arms.
- **Panels.** Nutrient and gut-hormone receptors, mechanosensors, vagal subtype markers, and
  fasting-switch genes with literature-predicted directions. *Calcr* and *Pmch* are absent from the
  deposited matrix.
- **Multiple testing.** Separate Benjamini–Hochberg families: A_nop, A_panel, A_ieg and C_panel.
- **Module B.** Spearman correlation of HFD-vs-fed and fasted-vs-fed log2 fold changes over genes
  with mean ≥ 10. The P value is exact over all 400 within-sex relabellings of HFD/fed animals; this
  null reproduces the positive correlation that shared controls induce. Descriptive split-control
  estimate: median ρ = 0.01 (range −0.35 to 0.37).
- **Early-life comparisons** reuse the PyDESeq2 fits from the primary analysis.

## Limits

- **Sample size:** 3 mice per sex per arm; bulk ganglion RNA mixes neuron subtypes with glia, so
  marker changes (*Uts2b*, *Vip*) may reflect cell composition rather than per-neuron regulation.
- **Normalized data only:** GSE281091 deposited normalized values only; raw counts from SRA would
  allow a count model with RNA-quality covariates.
- **No functional readouts:** these data cannot speak to food intake, vagal firing, energy
  expenditure or peptide levels.

## Figures (legends in `figures/LEGENDS.md`)

- **Fig. 4** Fasting response of the nodose ganglion.
- **Fig. 5** HFD vs fasting concordance and the satiety/nutrient-sensing panel across nodose datasets.
- **Fig. S4** Gene-length technical axis in GSE281091.
