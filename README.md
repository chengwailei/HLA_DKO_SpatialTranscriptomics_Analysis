# Immune and Mutational Profiling of Gene-Edited Low-Immunogenic Human Primary Cholangiocyte Organoids
This readme provides a brief description of replicating an example of the analysis included in the [Petrus-Reurer et al. ()] for the section on spatial transcriptomics analysis.

[![DOI]()

# Pacakge Installation
Ensure that R or RStudio is installed. The included example script requires three packages. 

The first two, Rtsne and umap are available from CRAN and can be installed with the following command in R:
``` r
install.packages(c("Rtsne", "umap"))
```
The following packages are managed by bioconductor and can be installed using the following command in R:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bio_pkgs <- c("NanoStringNCTools", "GeomxTools", "GeoMxWorkflows", "SpatialDecon", "standR", "limma", "ExperimentHub", "SpatialExperiment", "org.Hs.eg.db", "clusterProfiler")
BiocManager::install(bio_pkgs)
```
Typical installation time for all above packages: 45 mins

# Running the example
The included example R script reads the necessary intensity data out of the files in the `data` folder, producing the given output. To replicate the example you will need the following 5 files saved to your working directory:
- ./src/DKOvsWT_Analysis.Rmd
- ./src/DKO_ExpressionProfile.Rmd
- ./src/InVitro_Analysis.Rmd
- ./data/PlateA_GeoMx.rds
- ./data/PlateB_GeoMx.rds

## Bulk RNAseq Analysis
The analysis can be run from the RStudio directly from excuting the codes in the `InVitro_Analysis.Rmd` code, or using one of several options from R console. Once it is run, the following result contents will be generated:
- InVitro_Expression.csv
- InVitro_Sampledata.csv
- InVitro_ZScore_Dotplot.csv
- InVitro_Zscore_Expression.csv
  
## GeoMX Spatial Transcriptomics Analysis
The analysis can be run from the RStudio directly from excuting the codes in the `DKOvsWT_Analysis.Rmd` code, or using one of several options from R console. Once it is run, the following result contents will be generated:
- CD45_Deconvolution_SafeTME_Counts.csv: Deconvolution results of the CD45+ cells based on the SafeTME matrix.
- CD45_Deconvolution_ImmuneCensus_Counts.csv: Deconvolution results of the CD45+ cells based on the HCA Immune Census matrix.
- CD45_DEGs.csv: The differentially expressed genes between CD45+ cells within the compartment of DKO cholangiocytes and WT cholangiocytes.
- CD45_PathwayEnrichment.csv: The pathway enrichment results based on the differential expresssed genes between the CD45+ within compartment of DKO cholangiocytes and WT cholangiocytes.
- PanCK_DEGs.csv: The differentially expressed genes between the DKO cholangiocytes and WT cholangiocytes.
- PanCK_PathwayEnrichment.csv: The pathway enrichment results based on the differential expresssed genes the DKO cholangiocytes and WT cholangiocytes.

All of these results are also available in the `res` folder and can be used directly for the `Figure_Generation.Rmd` code to generate plots that are used in the paper.

Typical run time for analysis: 1 hour

# Directory structure
```
└─ src/:
|   └─ Figure_Generation.Rmd: Script to reproduce figures
|   └─ DKOvsWT_Analysis.Rmd: Analysis code used to reproduce results in paper
|   └─ DKO_ExpressionProfile.Rmd: Analysis comparing PanCK GeoMX (In Vivo e\PCOs) to primary cholangiocytes from Sampaziotis et al.
|   └─ InVitro_Analysis.Rmd: Analysis comparing in vitro e\PCOs expression.
└─ data/:
|   └─ PlateA_GeoMx.rds: Plate A data of GeoMx spatial transcirptomics
|   └─ PlateB_GeoMx.rds: Plate B data of GeoMx spatial transcirptomics
|   └─ InVitro_Expression.csv: log(TMM+1) expression of in vitro e\PCOs.
|   └─ InVitro_Sampledata.csv: In vitro e\PCOs sample information.
└─ res/:
|   └─ AllExpression_GeneralCOSummary.csv: General CO-associated expression of all In Vitro and In Vivo samples compared to priamry cholangiocytes from Sampaziotis et al.
|   └─ CD45_Deconvolution_SafeTME_Counts.csv: Deconvolution results of the CD45+ cells based on the SafeTME matrix.
|   └─ CD45_Deconvolution_ImmuneCensus_Counts.csv: Deconvolution results of the CD45+ cells based on the HCA Immune Census matrix.
|   └─ CD45_DEGs.csv: The differentially expressed genes between CD45+ cells within the compartment of DKO cholangiocytes and WT cholangiocytes.
|   └─ CD45_PathwayEnrichment.csv: The pathway enrichment results based on the differential expresssed genes between the CD45+ within compartment of DKO cholangiocytes and WT cholangiocytes.
|   └─ PanCK_DEGs.csv: The differentially expressed genes between the DKO cholangiocytes and WT cholangiocytes.
|   └─ PanCK_PathwayEnrichment.csv: The pathway enrichment results based on the differential expresssed genes the DKO cholangiocytes and WT cholangiocytes.
|   └─ InVitro_ZScore_Dotplot.csv: In Vitro e\PCOs expression by key regional genes.
|   └─ InVitro_Zscore_Expression.csv: In Vitro e\PCOs expression as normalized z-scores.
|   └─ InVivo_ExpressionHeatmap.csv: In Vivo e\PCOs (spatial transcriptomics) expression by key regional genes.
|   └─ InVivo_Prim_CorrelationMatrix.csv: Correlation expression of key regional genes between In Vivo e\PCOs (spatial transcriptomics) corerlation and primary cholangiocyte cells.
|   └─ InVivo_Prim_PCA_results.csv: Principal component analysis of key regional genes expressions between In Vivo e\PCOs (spatial transcriptomics) corerlation and primary cholangiocyte cells.
```

# Software Dependencies
The scripts were written using R (version 4.2.2) and the spatial transcriptomics analysis steps depend on the following pacakges:
- NanoStringNCTools (version 1.4.0)
- GeomxTools (version 3.0.1)
- GeoMxWorkflows (version 1.2.0)
- SpatialDecon (version 1.6.0)
- standR (version 1.8.0)
- limma (version 3.52.4)
- ExperimentHub (version 2.4.0)
- SpatialExperiment (version 1.6.1)
- Rtsne (version 0.17)
- umap (version 0.2.10.0)
- clusterProfiler (version 4.4.4)
- org.Hs.eg.db (version 3.15.0)
