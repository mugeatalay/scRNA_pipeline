##QC 

# List of dataset names and corresponding AnnData objects
adata_samples = {
    "control_1": adata1,
    "control_2": adata2,
    "control_3": adata3,
    "control_4": adata4,
    "treated_1": adata5,
    "treated_2": adata6,
    "treated_3": adata7,
    "treated_4": adata8,
}

for sample_name, adata in adata_samples.items():
    print(f"Generating QC metrics for {sample_name}...")

    # Ensure mitochondrial, ribosomal, and hemoglobin genes are labeled
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith("RPS") | adata.var_names.str.startswith("RPL")
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)  # Hemoglobin genes

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    # Generate and save QC plots
    sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    sc.pl.violin(adata, "pct_counts_mt", save=f"_violin_{sample_name}.pdf")
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

    print(f"Saved QC plots for {sample_name}.")

print("QC metrics and plots generated for all datasets.")

def preprocess_data(adata):
    # Ensure the AnnData is not a view
    adata = adata.copy()

    # Cell filtering
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
    mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = adata.X[:, mito_genes].sum(axis=1).A1 / adata.obs['n_counts'] * 100

    # Filter cells based on QC metrics
    adata = adata[(adata.obs['n_genes'] > 200) & 
                  (adata.obs['n_genes'] < 2500) & 
                  (adata.obs['n_counts'] < 50000) & 
                  (adata.obs['percent_mito'] < 5), :]

    # Filter ribosomal genes (RPS, RPL)
    ribosomal_genes_mask = ~(adata.var_names.str.upper().str.startswith('RPS') | 
                             adata.var_names.str.upper().str.startswith('RPL'))

    # Keep only non-ribosomal genes in the AnnData object
    adata = adata[:, ribosomal_genes_mask]

    return adata

# Preprocess individual datasets
    adata1 = preprocess_data(adata1) #control ZT0
    adata2 = preprocess_data(adata2) #control ZT0
    adata3 = preprocess_data(adata3) #control ZT12
    adata4 = preprocess_data(adata4) #control ZT12
    adata5 = preprocess_data(adata5) #treated ZT0
    adata6 = preprocess_data(adata6) #treated ZT0
    adata7 = preprocess_data(adata7) #treated ZT12
    adata8 = preprocess_data(adata8) #treated ZT12


## Check if any cell or gene contains infinite or NaN values
if np.any(np.isnan(adata1.X.toarray())) or np.any(np.isinf(adata1.X.toarray())):
    print("There are NaN or infinite values in the data. Please clean the data before proceeding.")
else:
    print("No NaN or infinite values detected.")

# Check for infinite values in the sparse matrix and remove them
if np.any(np.isinf(adata2.X.data)):
    print("Found infinite values in the data. Cleaning them...")
    # Replace infinite values with zeros or some other placeholder (like the max finite value)
    adata2.X.data[np.isinf(adata2.X.data)] = 0
else:
    print("No infinite values detected.")
