# Figure legends

Colours: control diet blue #0072B2 (circles); hypercaloric diet pink #CC79A7 (triangles).
Positive log2 fold changes indicate higher abundance under hypercaloric feeding. Confidence
intervals are nominal (not simultaneous). Source data: `source_data/`.

**Fig. 1 | Hypercaloric feeding does not reproducibly change *Pnoc* or *Oprl1* in peripheral ganglia.**
**a, b**, log2 fold change (hypercaloric/control) with nominal 95% CI for *Pnoc* (a) and *Oprl1* (b).
Each row is a prespecified within-study contrast; n C/H, independent animals (TRAP: pooled-DRG
immunoprecipitates) per arm. DRG (GSE225224, bulk; GSE198895, Nav1.8⁺ RiboTag TRAP) and nodose
(GSE320302) estimates: PyDESeq2 Wald tests on the full transcriptome. Nodose adult HFD (GSE281091):
moderated linear model on log2(DESeq2-normalized + 1) adjusted for sex, because raw counts are not
deposited. q, Benjamini–Hochberg across the primary family (2 genes × 5 contrasts with an estimate;
10 tests). Open squares, mean normalized count < 10. TRAP: both genes 0 counts in 4/4 libraries.
Sympathetic scRNA-seq (GSE233163): *Pnoc* absent from the deposited matrix; diet confounded with
sequencing run, so no diet effect is estimated. Rows are not pooled.

**Fig. 2 | Biological-replicate expression.** Normalized counts of *Pnoc* (a) and *Oprl1* (b) for each
animal; bar, median. Count series use PyDESeq2 median-of-ratios size factors from the full
transcriptome; GSE281091 shows the deposited DESeq2-normalized values (open symbols, females;
filled, males). No plotting transform. q as in Fig. 1. n as in Fig. 1. Note the different y scales.

**Fig. 3 | Cellular context for *Oprl1*.** **a**, Percentage of cells with ≥1 *Oprl1* UMI per deposited
cell type in sympathetic ganglia (GSE233163); one point per pool (site × sequencing run) with ≥20 cells.
**b**, *Oprl1* counts per million in neuron pseudobulk per pool (labels, sequencing run; lean pools from
N61/N66, obese pools from N68 — diet and run are confounded; descriptive only). Pools are of ganglia
from several mice. **c**, Percentage of putative neuronal (≥1 UMI *Snap25* and *Phox2b* or *Prdm12*)
and other nuclei with ≥1 UMI of *Pnoc* (open) or *Oprl1* (filled) in HFD-fed mouse vagus snRNA-seq
(GSE309905; Pten WT and vagal Pten KO; no control diet; libraries with ≥10 nuclei in the group).

**Fig. S1 | Sample quality control.** **a**, PCA of the 500 most variable genes (VST for count series;
log2(normalized + 1) for GSE281091); open symbols, later age; rings, median inter-sample correlation
robust z < −3.5 (flagged, retained). In GSE281091 PC1 isolates the three HFD-male samples and is
dominated by Gm/non-coding features. **b**, Library sizes (GSE281091: sum of deposited normalized values).

**Fig. S2 | Prohormone-processing enzymes.** As Fig. 1 for *Pcsk1*, *Pcsk2* and *Cpe* (separate BH
family; not specific readouts of N/OFQ signaling). Arrows, CI beyond axis.

**Fig. S3 | Exploratory opioid ligand/receptor panel.** As Fig. 1 for *Penk*, *Pdyn*, *Pomc*, *Oprm1*,
*Oprd1*, *Oprk1* (separate exploratory BH family). Numbers at the axis edge give estimates beyond the
axis; open squares, low abundance or detection in <50% of samples in an arm.

## Feeding-regulation figures (modules A–C; see `../REPORT_feeding_metabolism.md`)

Additional colours (Okabe–Ito): fasted 18 h orange #E69F00 (squares); fasted + refed 20 min bluish
green #009E73 (diamonds). Chow fed = control blue circles; HFD 4 wk = pink triangles.

**Fig. 4 | Fasting response of the mouse nodose ganglion (GSE281091).** **a**, Volcano plot,
18-h fasted vs ad libitum chow-fed (n = 6 vs 6, 3 M + 3 F each), moderated linear model on
log2(DESeq2-normalized + 1) adjusted for sex, 24-sample variance model; genes with mean normalized
count ≥ 10; coloured, transcriptome-wide BH padj < 0.05. **b–g**, Per-animal normalized counts for
*Pnoc*, *Oprl1* and prespecified feeding-panel genes across fed, fasted, refed (18-h fast + 20 min
food) and HFD (4 wk) mice; bar, median; open, females; ring, samples with length-bias index < −0.2
(technical axis, Fig. S4). q, BH within the prespecified family; length-adjusted P, post hoc
sensitivity model including the per-sample length-bias index.

**Fig. 5 | HFD does not reproduce the fasting state; robust HFD changes in vagal feeding genes.**
**a**, Gene-wise log2 fold changes, HFD vs fed against fasted vs fed (19,711 genes, mean ≥ 10;
hexagon density, log scale). **b**, Null distribution of Spearman ρ under all 400 within-sex
relabellings of HFD/fed animals (accounts for the shared fed controls); pink line, observed;
one-sided P. **c**, Satiety/nutrient-sensing panel: log2 fold changes for fasting (adult), HFD adult
(prespecified, female-only and length-adjusted sensitivities) and early-life HFD at P21/P35
(GSE320302, PyDESeq2, n = 8/8 and 8/7). Filled dot, BH q < 0.05 within the prespecified family;
open dot, nominal P < 0.05; ×, absent from matrix or mean < 5. **d–f**, Per-animal values for the
three HFD effects robust to both sensitivities. Symbols as Fig. 4.

**Fig. S4 | Technical axis in GSE281091.** **a**, PCA of the 500 most variable expressed genes
(24 samples); PC2 separates sex. **b**, PC1 against the per-sample length-bias index (Spearman ρ
between centred log expression and log10 gene length across expressed genes). The three HFD-male
samples carry strong length bias, typical of reduced RNA integrity or 3′ bias.

**Fig. 6 | No feeding-regulatory gene is consistently changed by HFD across ganglia.**
**a**, HFD log2 fold changes in DRG (GSE225224; inverse-variance combination of 12- and 36-wk
PyDESeq2 estimates) vs nodose (GSE281091; moderated model adjusted for sex and length-bias index)
for 184 expressed feeding-regulatory genes (GO feeding/eating/appetite/energy-homeostasis/leptin
sets plus curated vagal panel); orange, fasting (nodose, 18 h) changes in the opposite direction to
nodose HFD; grey dots, all other genes. No gene reached BH q < 0.1 in the conjunction test
(changed in both ganglia). **b**, The 16 genes with the smallest conjunction P plus *Lepr*, *Uts2b*,
*Pnoc*, *Oprl1*, across all datasets: DRG 12 and 36 wk (n = 6/6, 5/6), Nav1.8⁺ DRG TRAP (2/2 pools),
nodose adult HFD (6/6) and early-life HFD at P21 (8/8), sympathetic SCG and stellate obese pool vs
lean pools (descriptive; diet confounded with sequencing run; dot, obese–lean difference no larger
than the lean between-run difference), and nodose fasting. Row label: consensus direction and
conjunction P (– if P = 1). ○, nominal P < 0.05. **c**, Adipocyte markers in DRG and nodose HFD
contrasts; the DRG 12-wk *Pparg* increase is accompanied by *Cidea*/*Cidec*/*Fabp4*, suggesting
more perineural adipose tissue in the dissected samples.
