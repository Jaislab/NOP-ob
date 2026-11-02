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

---

# Which feeding-regulatory genes does HFD change in all ganglia, ideally opposite to fasting?

Prespecified in `metadata/hfd_consensus_plan.md` (commit 86ebb21). Results:
`results/consensus/feeding_gene_consensus.csv`; Fig. 6.

## Answer

**No feeding-regulatory gene is reproducibly changed by HFD across ganglia, and none shows
significant HFD/fasting counter-regulation.**

- **Gene universe:** 289 genes from GO feeding/eating/appetite/energy-homeostasis/leptin sets plus
  the curated vagal panel; 184 are expressed in both primary ganglia (DRG, GSE225224; nodose,
  GSE281091).
- **Conjunction test** ("changed in DRG *and* nodose"): best P = 0.024, every q = 1.0.
  Transcriptome-wide, 0 of 13,947 genes reach q < 0.1.
- **Detectable effects:** the design could detect shared changes of ≤ 0.42 log2 (≈ 34%)
  in half of the feeding genes and ≤ 0.25 log2 in a quarter (80% power, one-sided α = 0.05 in each ganglion). Shared
  regulation larger than that is unlikely for most well-expressed feeding genes; smaller or
  cell-type-restricted changes cannot be excluded.

## Best nominal candidates (not significant after correction)

| Gene | HFD direction | DRG 12/36 wk | Nodose HFD | Fasting (nodose) | Comment |
|---|---|---|---|---|---|
| *Pparg* | ↑ | +1.60 / +0.60 | +0.49 | +0.07 | Also ↑ in Nav1.8 TRAP and early-life nodose; in DRG it rises with adipocyte genes (*Cidea*, *Cidec*, *Fabp4*), so extra fat around the ganglia is a likely contributor; not counter-regulated by fasting |
| *Ugcg* | ↑ | +0.07 / +0.24 | +0.15 | −0.01 | Small effects |
| *Acacb* | ↓ | −0.33 / −0.31 | −0.23 | **+0.27** | Counter-regulated by fasting (P_counter = 0.12) |
| *Insr* | ↓ | −0.15 / −0.09 | −0.12 | **+0.09** | Counter-regulated (P_counter = 0.14) |
| *Cnr1* | ↓ | −0.09 / −0.14 | −0.11 | **+0.21** | Counter-regulated (P_counter = 0.25) |
| *Pomc* | ↑ | +0.79 / +0.30 | +0.19 | **−0.26** | Counter-regulated (P_counter = 0.34) |
| *Cartpt* | ↑ | +0.60 / −0.32 | +0.23 | **−0.11** | Counter-regulated (P_counter = 0.39) |

- **A pattern worth noting:** the counter-regulated candidates make a coherent "energy-surplus"
  pattern. Anorexigenic *Pomc*/*Cartpt* go up with HFD and down with fasting; orexigenic/anabolic
  *Cnr1*, *Insr* and *Acacb* go down with HFD and up with fasting. The effects are small (≤ 0.3
  log2 in nodose) and individually not significant.
- **The pattern is not enriched.** Among feeding genes changing in the same direction in DRG and
  nodose, 45% are counter-regulated by fasting, the same as transcriptome-wide (45%; Fisher
  P = 1.0). It may be chance.

## Ganglion-specific rather than shared

- **Nodose HFD changes stay in the nodose.** The robust nodose changes (*Uts2b* ↑, *Vip* ↑,
  *Lepr* ↓) do not appear in DRG: *Lepr* +0.24, *Uts2b* −0.23, and *Vip* is not expressed in DRG.
- **DRG *Pnoc* ↓ stays in the DRG.** The 12-wk DRG *Pnoc* decrease is absent in nodose (−0.05).

## Sympathetic ganglia (GSE233163, descriptive)

- **Confounded design:** whole-ganglion pseudobulk from the obese pool (run N68) was compared with
  lean pools (runs N61, N66).
- **Run noise matches the diet differences:** for most candidates the obese–lean difference is no
  larger than the difference between the two lean runs (*Pparg* SCG −0.94 vs a lean run difference
  of −0.78), and *Pparg* goes the opposite way to DRG/nodose.
- **No conclusion is possible**, because diet is confounded with run.

## Caveats

- **Fasting data:** fasting data exist for nodose only. GEO searches (2 Nov 2026) found no DRG or
  sympathetic-ganglion fasting datasets: GSE239596 is an insulin-receptor-antagonist model and
  GSE224316 is DRG leptin overexpression.
- **Nodose model choice:** the primary nodose model adjusts for the RNA-quality (gene-length) axis.
  Prespecified-model and female-only versions of every test are in the results table.
- **Different exposures:** DRG (45% HFD, 9–33 wk, males) and nodose (60% HFD, 4 wk, both sexes)
  differ in diet, duration and sex, which reduces the chance of seeing shared effects.

---

# Are the HFD/fasting counter-regulated candidates expressed in the geniculate ganglion?

Fig. 7; `results/geniculate/`. Expression only: no diet or fasting data exist for the geniculate
ganglion in GEO (searched 2 Nov 2026).

**Datasets**

- **GSE223724:** bulk geniculate ganglia, 4 wild-type adults; the *Egr4*⁻/⁻ mice were not used.
- **GSE102443:** 96 single geniculate neurons (Fluidigm C1, FPKM), split by *Phox2b*:
  - 61 *Phox2b*⁺ oral/taste neurons
  - 35 *Phox2b*⁻ somatosensory neurons
- **GSE135801:** 376 *Phox2b*⁺ geniculate neurons with ≥ 100,000 reads; this cutoff was fixed
  before looking at the target genes.
- **References:** nodose (chow fed) and DRG (standard diet), converted to TPM with the same gene
  lengths.

| Gene | Bulk geniculate TPM (4 WT) | Nodose / DRG TPM | Neurons expressing: C1 *Phox2b*⁺ / *Phox2b*⁻ / GSE135801 | Verdict |
|---|---|---|---|---|
| *Acacb* | 1.1 (0.5–1.7) | 4.8 / 3.4 | 3% / 0% / 5% | Low; essentially not neuronal (likely glial/vascular) |
| *Insr* | 1.7 (1.0–2.4) | 17.5 / 9.8 | 72% / 71% / 31% | Low in bulk, but widely detected at low level in neurons |
| *Cnr1* | 6.5 (5.0–7.4) | 94.6 / 19.9 | **97%** / 46% / 80% | Expressed, mainly in oral/taste (*Phox2b*⁺) neurons |
| *Pomc* | 2.6 in the one pituitary-free sample* | 2.6 / 1.7 | 43% / 31% / 4% (≈ 2 FPKM) | Trace only |
| *Cartpt* | 17.0 (11.6–21.5) | 368 / 63.5 | 51% / **100%** / 23% | Expressed; in all somatosensory (*Phox2b*⁻) neurons and about half of the taste neurons |

\*Bulk geniculate *Pomc* ranged from 2.6 to 1,381 TPM and tracked the pituitary hormones *Gh*,
*Prl* and *Cga* sample by sample (pituitary carry-over from the skull-base dissection). Only one
wild-type sample was free of pituitary signal (`geniculate_pituitary_check.csv`).

**Summary**

- **Expressed:** *Cnr1* and *Cartpt* are clearly expressed by geniculate neurons, with the same
  split as in the vagal ganglion: CB1 in taste/visceral-type (*Phox2b*⁺) neurons, and CART strongest
  in somatosensory neurons.
- **Low:** *Insr* is present at a low level in many neurons.
- **Effectively absent:** *Acacb* and *Pomc*.
- **Lower bulk levels:** both genes are lower in the geniculate than in the nodose (about 15-fold for
  *Cnr1* and 20-fold for *Cartpt*). The geniculate is a small, glia-rich ganglion; *Snap25* is
  similar across ganglia, so part of the difference is biological, but cross-study TPM comparisons
  are approximate.
- **If you go further:** the geniculate is a reasonable place to test *Cnr1* and *Cartpt* under HFD
  or fasting. *Pomc*, *Acacb* and probably *Insr* are poor geniculate readouts.
