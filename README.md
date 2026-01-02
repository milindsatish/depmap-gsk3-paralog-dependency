# depmap-gsk3-paralog-dependency
DepMap-based analysis of GSK3A and GSK3B paralog dependencies across cancer cell lines. Integrates CRISPR gene effect data with transcriptomic profiles to identify context-specific, paralog-selective dependencies and associated biological programs relevant to target prioritization.

# Context-Specific and Paralog-Selective Dependencies in DepMap: A GSK3 Case Study

## Overview

• Identify transcriptional programs associated with GSK3B dependency.
• Lay groundwork for precision biomarker discovery.

This project demonstrates how integrated functional genomics data can be used to identify **context-specific, paralog-selective target dependencies** and link them to underlying transcriptional programs. Using the GSK3 paralogs (GSK3A and GSK3B) as an illustrative example, the analysis highlights how gene dependency patterns vary across cellular contexts and how these differences can inform **rational target prioritization and patient stratification**.

The analysis leverages data from the Cancer Dependency Map (DepMap) to examine whether gene expression patterns are associated with dependency on the kinase **GSK3B** across human cancer cell lines. Gene expression, CRISPR gene dependency, and sample metadata are integrated into a single aligned dataset, with an emphasis on **data integrity, reproducibility, and biological interpretability** rather than black-box prediction.

Work completed to date establishes a rigorous multi-omics data integration pipeline and performs exploratory, hypothesis-generating analyses that serve as a foundation for downstream modeling and translational interpretation.

---

## Data Sources

Three publicly available DepMap datasets are used:

1. **CRISPR Gene Effect (Chronos)**  
   Genome-wide CRISPR knockout dependency scores, where more negative values indicate stronger gene essentiality.

2. **CCLE Gene Expression (Affymetrix microarray)**  
   Probe-level gene expression data measured across cancer cell lines.

3. **Sample Metadata (`sample_info.csv`)**  
   Cell line annotations including DepMap identifiers, tissue lineage, and disease subtype.

---

## Data Loading and Initial Quality Control

All datasets are loaded using `data.table` for efficiency and inspected to confirm dimensions, structure, and data integrity. Expression data are read from GCT format after skipping header lines, and initial previews verify gene annotation columns, sample identifiers, and numeric consistency.

---

## Mapping Expression Samples to DepMap Identifiers

CCLE expression sample names are mapped to DepMap cell line identifiers using a lookup table derived from sample metadata. Expression columns are renamed from CCLE names to DepMap IDs, duplicated samples are removed, and multiple validation checks confirm identifier completeness and uniqueness. After mapping, **988 expression samples** are retained.

---

## Sample Intersection and Alignment

The three datasets differ in orientation and identifier placement: expression data store samples in columns, CRISPR data store samples in rows, and metadata store samples in rows. The intersection across all datasets yields **677 shared cell lines**. Expression is the limiting dataset; all retained expression samples have corresponding CRISPR dependency data and metadata.

All datasets are subset and reordered so that the same sample appears in the same position throughout. Explicit sanity checks confirm identical sample ordering across expression, CRISPR, and metadata tables, eliminating the possibility of silent misalignment.

---

## Target Definition: GSK3B Dependency

The biological target of interest is **GSK3B (Glycogen Synthase Kinase 3 Beta)**. The CRISPR dependency score for GSK3B is extracted as the dependent variable, yielding one numeric value per cell line that quantifies the effect of GSK3B knockout on cell viability. The resulting distribution shows substantial heterogeneity across cell lines, indicating **context-specific essentiality** and motivating downstream analyses.

---

## Feature Matrix Construction

Gene expression values are used as independent variables. Gene annotation columns are removed, the expression matrix is transposed so that rows correspond to samples and columns correspond to genes, and gene identifiers are explicitly reassigned after transposition. Alignment between expression features and the dependency vector is verified to ensure sample-level consistency.

---

## Exploratory Correlation Analysis

As an interpretable, hypothesis-generating first step, **Pearson correlations** are computed between expression of each gene and GSK3B dependency across the 677 aligned cell lines. Correlations are calculated gene-by-gene, producing a ranked list of association strengths.

Probe IDs are mapped to gene symbols using expression dataset annotations, and Affymetrix control probes are removed. The resulting correlation distribution is centered near zero, with most genes showing minimal association and a subset exhibiting modest positive or negative correlations (approximately −0.2 to +0.2). This pattern is consistent with **distributed, context-dependent regulation** rather than a single dominant transcriptional driver.

All results at this stage are exploratory and associative rather than causal.

---

## GSK3A vs GSK3B Dependency Comparison

To assess whether GSK3B dependency reflects general GSK3 biology or paralog-specific effects, CRISPR dependency scores for **GSK3A** and **GSK3B** are compared directly across all aligned cell lines. For each cell line, GSK3A dependency is plotted against GSK3B dependency, and the relationship is quantified using Pearson correlation.

Across all cell lines, dependency on GSK3A and GSK3B is **moderately correlated** (Pearson r ≈ 0.34), indicating shared but incomplete functional overlap between the two paralogs. Stratification by tissue lineage shows substantial overlap across lineages, while off-diagonal dispersion indicates **paralog-selective and context-dependent essentiality** rather than simple redundancy.

---

## Identification of GSK3B-Selective Cell Lines

To identify paralog-selective vulnerabilities, a **GSK3B selectivity score** is defined as the difference between GSK3B and GSK3A dependency scores for each cell line. While most cell lines show minimal selectivity, a subset exhibits strong GSK3B-specific dependency.

Ranking cell lines by this metric identifies **68 GSK3B-selective cell lines** (bottom 10%), spanning multiple tissue lineages. This supports the presence of context-dependent, paralog-specific dependencies and motivates pathway-level analyses.

---

## Pathway-Level Enrichment in GSK3B-Selective Cell Lines

Pathway enrichment analysis is performed to identify transcriptional programs associated with GSK3B-selective dependency using ranked gene-level associations.

**Top enriched pathways** include hormone response and metabolic programs, such as:
- Estrogen Response (Early and Late)
- Cholesterol Homeostasis
- Interferon Alpha Response
- Androgen Response
- Notch Signaling

In contrast, **depleted pathways** include canonical proliferation-associated programs, such as:
- WNT / β-catenin signaling
- MYC targets
- E2F targets
- G2/M checkpoint
- Mitotic spindle

Together, these results indicate that **GSK3B-selective dependency is associated with hormone-responsive and differentiation-related transcriptional states**, rather than generalized proliferative or WNT-driven essentiality. This distinction is directly relevant for **target prioritization and potential patient stratification** in drug discovery.

---

## Validation of Estrogen Response Pathway Activity

To validate pathway enrichment results, estrogen response pathway activity is quantified per cell line as the mean expression of Hallmark estrogen response (early and late) genes. GSK3B-selective cell lines show a higher median estrogen response score compared with other cell lines, with a consistent shift in the distribution rather than a small number of extreme outliers. This supports a robust association between **GSK3B-selective dependency and hormone-responsive transcriptional states**.

---

## Summary and Outlook

This project illustrates how integrated analysis of DepMap CRISPR dependency and expression data can uncover **context-specific and paralog-selective vulnerabilities** with clear biological interpretation. The results highlight GSK3B-selective dependency as a context-dependent phenotype associated with hormone response and metabolic programs, providing a framework for downstream modeling, target refinement, and translational hypothesis generation.
