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
