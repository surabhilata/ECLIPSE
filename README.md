# ECLIPSE
The ECLIPSE (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) framework embeds organism-specific proteomes into a global protein sequence similarity network derived from the Protein Universe Atlas
Repo organisation

# Repo organisation
The code is organised in python notebooks (for dark proteome data analysis), python scripts for (cluster analysis). ECLIPSE contains two notebooks each notebook is devided in two major analysis. 

 # Dependencies
The code was written in Python 3.12.4.

For the analysis of the data the modules required are:

    1. jupyter
    2. pandas
    3. seaborn
    4. requests
    5. biophython
    6. jedi==0.17.2
    7. mongodict
    8. memory_profiler
    9. pymongo
    10. holoviews
    11. datashader
    12. holoviews[recommended]
    13. umap-learn
    14. faerun

# How to use ECLIPSE

# Initiate the notebook environment by command ("mamba activate venv")

# To perform large scale sequence similarity search

1. We have used MMseq2 easy-search where we have used P.aeruginosa panproteome dataset in tar file which was searched over Atlas AFDBv4_90.fasta file. #For each protein input file we have taken the best match with following command in terminal
    # mmseqs easy-search PA.tar.gz AFDBv4_90.fasta atlas_search_results.m8 tmp --max-seqs

2. 
      

     

 
