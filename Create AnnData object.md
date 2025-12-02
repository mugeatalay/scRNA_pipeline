# Import Modules 
## 1.1 Install Libraries

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
import seaborn as sns
import leidenalg
import scipy.sparse as sparse
import scipy.io
import scrublet as scr
import session_info

## 1.1.1 Load Data
## 1.1.2 Read Data

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi = 160, color_map = 'RdPu', dpi_save = 180, vector_friendly = True, format = 'svg')
session_info.show()


Sesh info**** put

def charge_adata(matrix_file, barcode_file, feature_file):

    # Load the count matrix
    counts_matrix = scipy.io.mmread(matrix_file).T.tocsc()  # Transpose to have cells as rows and genes as columns
    barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')
    features = pd.read_csv(features_file, header=None, sep='\t')
    
    # Create AnnData object
    adata = sc.AnnData(X=counts_matrix)
    adata.obs['barcode'] = barcodes[0].values
    adata.var['gene_names'] = features.iloc[:, 1].values  # Select only the second column (gene names)
    adata.var['ensembl_id'] = features.iloc[:, 0].values  # Select the Ensembl IDs (1st column) if needed later
    
    # Set gene names as variable names and ensure they are unique
    gene_names = features.iloc[:, 1].values  # Assuming the second column contains gene names
    
    # Assign gene names and make them unique to avoid issues with duplicated names
    adata.var_names = pd.Index(gene_names).astype(str)
    adata.var_names_make_unique()

    return adata

# Load your data as AnnData (see previous examples for loading data)
matrix_file = '/data/scc/LIPIMMUNE_05/A-PCB-24014554_L_GEX_folder/matrix.mtx' # Your count matrix
barcodes_file = '/data/scc/LIPIMMUNE_05/A-PCB-24014554_L_GEX_folder/barcodes.tsv'  # Your cell barcodes
features_file = '/data/scc/LIPIMMUNE_05/A-PCB-24014554_L_GEX_folder/features.tsv'  # Your gene features

adata1 = charge_adata(matrix_file, barcodes_file, features_file)

Repeat this as number of samples you wish to analyze
