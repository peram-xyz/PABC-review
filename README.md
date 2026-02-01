 # PABC-review

Data File Sources:
 - Transcriptomic profiles pregnancy-associated breast cancer (PABC) tumours were downloaded from the GEO database (https://www.ncbi.nlm.nih.gov/geo/) from a study characterizing the molecular signature of PABC under accession code GSE31192
 - Transcriptomic profiles of breast invasive carcinoma tumours were downloaded from the TCGA-BRCA project (https://portal.gdc.cancer.gov/projects/TCGA-BRCA)
 - Transcriptomic profiles of placental cells were also downloaded from the GEO database from a study conducting single-cell RNA-seq analysis on trophoblast subtypes under accession code GSE89497

Data files are omitted from this repository, but can be downloaded from the respective sites above.

The PyOD, pyCombat, and mygene libraries were used for data analysis. Furthermore, the pyCirclize library was used to visualize pathway enrichment analysis results.

![chord plot](data/processed/circos_plot_1003.png)

All analyses were run on Python (v 3.13.5).
