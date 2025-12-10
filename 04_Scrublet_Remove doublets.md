
## Import
```python
import scrublet as scr

def detect_doublets(adata):
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=30)
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets
    
    # Filter out doublets
    adata = adata[~adata.obs['predicted_doublets'], :]
    return adata

# Run Scrublet on multiple datasets
adata1 = detect_doublets(adata1)
adata2 = detect_doublets(adata2)
adata3 = detect_doublets(adata3)
adata4 = detect_doublets(adata4)
adata5 = detect_doublets(adata5)
adata6 = detect_doublets(adata6)
adata7 = detect_doublets(adata7)
adata8 = detect_doublets(adata8)
```

<img width="432" height="187" alt="DoubletScore_im" src="https://github.com/user-attachments/assets/4879d8ef-5dd5-450e-83ad-c4ff453f7bc0" />

## Plot Scrublet doublet scores
```python
import matplotlib.pyplot as plt

def plot_doublets(adata, title="Doublet Score Distribution"):
    plt.figure(figsize=(6,4))
    plt.hist(adata.obs['doublet_score'], bins=50)
    plt.xlabel("Doublet Score")
    plt.ylabel("Number of Cells")
    plt.title(title)
    plt.show()


# Example: plot for each dataset
plot_doublets(adata1, "Dataset 1 Doublet Scores")
plot_doublets(adata2, "Dataset 2 Doublet Scores")
```
scrub.plot_histogram()

Hist*





