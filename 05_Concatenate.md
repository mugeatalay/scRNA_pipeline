## **Create condition labels**

```python
list1 = ['Young 1'] * len(adata1)  # batch no 3
list2 = ['Young 2'] * len(adata2)  # batch no 1
list3 = ['Young 3'] * len(adata3)  # batch no 3
list4 = ['Young 4'] * len(adata4)  # batch no 2
list5 = ['Aged 1'] * len(adata5)   # batch no 1
list6 = ['Aged 2'] * len(adata6)   # batch no 2
list7 = ['Aged 3'] * len(adata7)   # batch no 2
list8 = ['Aged 4'] * len(adata8)   # batch no 1
```

# Concatenate condition labels
label = np.concatenate([list1, list2, list3, list4, list5, list6, list7, list8])

```python
# Define the real batch numbers corresponding to each dataset
batch_numbers = (
    [3] * len(adata1) +  # A-PCB-24014554 
    [1] * len(adata2) +  # 240009245
    [3] * len(adata3) +  # A-PCB-24010581
    [2] * len(adata4) +  # 240009242
    [1] * len(adata5) +  # 220049299
    [2] * len(adata6) +  # 220049300
    [2] * len(adata7) +  # 220049298
    [1] * len(adata8)    # 220049297
)
```

# Concatenate AnnData objects with a neutral batch key (not "batch")
adata = adata1.concatenate(
    adata2, adata3, adata4, adata5, adata6, adata7, adata8, 
    batch_key='dataset_id',  # Avoid naming it 'batch'
    join='outer'
).copy()

# Assign the condition labels
adata.obs['sample'] = label  

# Assign real batch numbers
adata.obs['real_batch'] = batch_numbers  # This will be used for Harmony

# Create replicate labels
replicate_labels = (
    ['replicate1'] * len(adata1) +
    ['replicate2'] * len(adata2) +
    ['replicate1'] * len(adata3) +
    ['replicate2'] * len(adata4) +
    ['replicate1'] * len(adata5) +
    ['replicate2'] * len(adata6) +
    ['replicate1'] * len(adata7) +
    ['replicate2'] * len(adata8)
)

# Add the 'replicate' column
adata.obs['replicate'] = replicate_labels  
