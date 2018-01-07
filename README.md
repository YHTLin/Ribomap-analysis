# Ribomap-analysis
Analysis of ribosome profiling data on the effect of mTOR inhibitors in GTML5 mouse cells.

## Motivation
Signaling pathways are responsible for coordinating cell behavior, and when dysregulated they often lead to
disease. The PIK3CA-AKT-mTOR pathway is commonly activated in human cancers. It plays an important
role in cell cycle regulation, metabolism, and protein synthesis. One of its members, the mechanistic target
of rapamycin (mTOR), is an serine/threonine kinase whose over-activation promotes tumor progression.
To reverse the effect of increased mTOR activity, mTOR inhibitors, such as rapamycin (first generation)
and mTOR kinase inhibitors (second generation), are prescribed as cancer therapy. However, over time,
natural selection of tumors results in drug resistance and clinical relapse. To counteract mTOR mutants,
Rodrik-Outmezguine et al. recently developed a third-generation mTOR inhibitor called Rapalink-1 (M1071).
M1071 has been established by Fan et al. as a potent drug for brain tumors owing to its ability to cross the
blood-brain barrier. While the mechanism of action of M1071 is known by design, its downstream effects on
biological processes, such as protein translation, are not well understood. The goal of this study is to identify
proteins whose biosyntheses are sensitive to M1071 in the context of neural cells. We aim to monitor changes
in protein translation using two techniques, ribosome profiling and pulsed SILAC proteomics.

## Instructions
1. Download the **STATS** file containing the Ribomap output from the [Kingsford Group](https://github.com/Kingsford-Group/ribomap).
2. Download the R scripts and update the directory in the script to the location of the STATS file.
3. Run **ribosome_profiling_analysis_technical_rep.R** followed by **ribosome_profiling_report_ozlem.R** to generate a comprehensive report.

*NOTE: **ribosome_profiling_analysis_biological_rep.R** performs the statistical tests assuming that the replicates are biologically related (paired). However, the replicates in this experiment were gathered on the same day. The biological_rep script has since been updated and saved as the technical_rep version, which addresses this subtlety and implements an additional filtering function.*
