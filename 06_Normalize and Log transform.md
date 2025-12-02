```python
# Make an explicit copy of your AnnData object to avoid warnings 
adata.layers['counts'] = adata.X.copy()

# Normalize the data (total-count normalization and log transformation)
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize
sc.pp.log1p(adata)  # Log-transform data

# Save AnnData after normalization
adata.write('Anndata0.h5ad')
```
