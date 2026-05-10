# Paper title

Co-expression network analysis reveals repression of Igf2bp2 by REV-ERBβ in skeletal muscle

## Contents

- Python-based network and correlation analysis
- R-based RNA-seq and WGCNA analysis
- Circadian dataset visualization in python
- Python notebook for generating some of the figures
- R-based code for ChIP-seq data plotting


## Folder structure

├── LICENSE
├── python_environment.yml
├── Python_scripts
│   ├── 01_Network_analysis
│   │   └── network_analysis_turquoise.ipynb
│   ├── 02_figure_1_and_2_panels
│   │   └── figure_1_and_2_panels.ipynb
│   └── 03_circadian_datasets_plotting_and_correlation
│       └── circadian_datasets_plotting_expression_and_correlation.ipynb
├── README.md
├── r_environment.yml
└── R_scripts
    ├── 01_rnaseq_preprocessing
    │   ├── RNAseq_preprocessing.md
    │   └── salmon_quants_to_tximport_wgcna.R
    ├── 02_wgcna_analysis
    │   └── WGCNA_GSE197726.R
    ├── 03_Differential_gene_expression
    │   └── Differential_gene_expression.R
    ├── 04_Differential_rhythmicity_and_rhythmicity_analysis
    │   ├── Differential_rhythmicity_DryR.R
    │   ├── metacycle_analysis_circadian_datasets.R
    │   ├── rain_analysis_circadian_datasets.R
    │   └── Rain_and_metacycle_output_heatmap.R
    ├── 05_chipseq_data_plotting
    │   └── chip_Seq_data_plotting.R
    └── 06_common_neighbors_and_circadian_community_GO
        └── common_neighbors_and_circadian_community_GO.R

## Authors and contact information

### Authors
A.M., Vishnu1†; Zhang, Qing2†; Srivastava, Shriyansh1; Koronowski, Kevin B.2*; Srivastava, Ashutosh1*

1Department of Biological Sciences and Engineering, Indian Institute of Technology Gandhinagar, Gandhinagar, Gujarat, India

2Department of Biochemistry & Structural Biology, University of Texas Health San Antonio, San Antonio, Texas, USA
Sam and Ann Barshop Institute for Longevity and Aging Studies, University of Texas Health San Antonio, San Antonio, Texas, USA

†Equal contribution

*Co-corresponding authors. 

### Lead contact
Dr. Ashutosh Srivastava (codes) and Dr. Kevin B. Koronowski (resources and reagents).

ashutosh.s@iitgn.ac.in, koronowski@uthscsa.edu 

