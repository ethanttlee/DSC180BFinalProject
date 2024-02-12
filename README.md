### Coloc Matrix Generation

We generate the coloc matrix in R using the code containted in the "coloc_matrix.ipynb" file iterating through disease datasets and cis-eQTL summary statistics generated with 1000 genomes. Reproducing the resulting matrix requires the following:

1. Generate cis-eQTL summary statistics by running "generate_eQTL_base_data.py", with the 1000 G European Individuals downloadable at https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 , gene expression data for the individuals "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt" and gene annotation file "gene_annot.txt" in the same directory as the script.
2. Downloading GWAS data using the "load_gwas.sh" bash script
3. Running "coloc_matrix.ipynb". This step requires that you install the R package coloc.

Since these operations are all computationally intensive, it is recommended that they be performed in a remote development server. A docker environment suitable for all tasks with coloc installed can be found at ghcr.io/jjdrisco/dsmlp-coloc-notebook .


### Gene Set Enrichment Analysis


The first part of our project is the gene set enrichment analysis, which we performed on gene expression data from the 1000Genomes dataset and the Hallmark gene sets. Our analysis began with the data manipulation, as seen in the "GSEA.ipynb" Jupyter notebook.

- Download the entire repository (including the data files), and then open the notebook and run it.
- Follow each cell step by step, and by the end, the appropriately formatted data should be located in a new file called "gene_exp_format_common.txt".
- Download this file locally, as it will be an input into the GSEA software.

- Install the GSEA software here: https://www.gsea-msigdb.org/gsea/downloads.jsp.
- Follow download instructions, and open the application.
- Under the 'Load data' tab, upload the gene expression data ("gene_exp_format_common.txt") and the label data ("age.cls").
    - The gene expression data should be generated by the notebook locally, and age.cls is included in the github.
- Go to the 'Run GSEA' tab, and load the following data
    - expression dataset: "gene_exp_format_common.txt"
    - phenotype labels: "age.cls"
- Select the gene sets database
    - Hallmark gene set: ("h.all.v2023.2.Hs.symbols.gmt")
- Select the menu for the "Collapse/Remap to gene symbols"
    - change value to "No_Collapse"

If the GSEA is completed, you will see "Running" then "Success" ath the bottom left of the screen.
These results include enrichment scores and charts for each gene set, as well as data for each individual gene. 
