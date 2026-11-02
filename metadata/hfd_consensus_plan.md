# Cross-ganglion HFD consensus of feeding-regulatory genes, with fasting counter-regulation

Frozen before estimation (2 Nov 2026). Question: which feeding-regulatory genes change with
HFD feeding in all ganglia, preferably in the direction opposite to fasting?

## Gene universe (feeding-regulatory)
Union of MSigDB 2024.1.Mm GO:BP sets (m5.go.bp.v2024.1.Mm.symbols.gmt, sha256 07e91b5e...):
FEEDING_BEHAVIOR, ADULT_FEEDING_BEHAVIOR, EATING_BEHAVIOR, REGULATION_OF_FEEDING_BEHAVIOR,
POSITIVE/NEGATIVE_REGULATION_OF_FEEDING_BEHAVIOR, REGULATION_OF_EATING_BEHAVIOR,
POSITIVE/NEGATIVE_REGULATION_OF_EATING_BEHAVIOR, REGULATION_OF_APPETITE,
REDUCTION_OF_FOOD_INTAKE_IN_RESPONSE_TO_DIETARY_EXCESS, ENERGY_HOMEOSTASIS, RESPONSE_TO_LEPTIN,
CELLULAR_RESPONSE_TO_LEPTIN_STIMULUS, LEPTIN_MEDIATED_SIGNALING_PATHWAY (all GOBP_), plus the
curated panels in feeding_gene_panels.csv (vagal_feeding, nop). Testable = mean normalized count
>= 10 in every primary dataset.

## Datasets
Primary (adult, ongoing HFD, within-study controls):
1. DRG, GSE225224: HFD vs SD at 12 wk and 36 wk (PyDESeq2 Wald), combined by inverse-variance
   fixed effect into one DRG estimate.
2. Nodose, GSE281091: HFD vs fed chow, moderated model adjusted for the per-sample length-bias index
   (the technical axis found earlier). Sensitivities: prespecified 24-sample model; female-only.
Supportive (reported, not required for "all"): Nav1.8+ DRG TRAP GSE198895 (n = 2/arm); early-life
nodose GSE320302 P21 (developmental exposure); sympathetic SCG/stellate GSE233163 (descriptive only:
whole-ganglion pool pseudobulk, obese pool vs lean pools per site, compared with the lean
between-run difference; diet confounded with run).

## Tests
- HFD in all primary ganglia: intersection-union (conjunction) test. For each direction,
  p_dir = max of the one-sided P values in DRG and nodose; p_conj = min(1, 2 x min over directions).
- HFD consensus counter-regulated by fasting: conjunction adding nodose fasted vs fed
  (length-adjusted) one-sided in the opposite direction: p_counter = min(1, 2 x min over directions
  of max(p_DRG,dir; p_nod,dir; p_fast,opposite)).
- BH across testable universe genes, separately for p_conj and p_counter. Report q < 0.1 and
  rank lists. Shared fed controls make HFD and fasting estimates positively correlated, so the
  opposite-sign requirement is conservative.
