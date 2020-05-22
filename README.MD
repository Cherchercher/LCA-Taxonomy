# LCA Taxonomy

This python program finds the nearest common ancestors of two specie name

### Setup

get data set to current work directory
`wget "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
`tar -xvf taxdump.tar.gz`

install packages
`pip install requirements.txt`

### Usage:

python lca.py name1 name2

### Example Execuation

```
python lca.py Monera Bacteria
reading files....
lca.py:125: ParserWarning: Falling back to the 'python' engine because the 'c' engine does not support regex separators (separators > 1 char and different from '\s+' are interpreted as regex); you can avoid this warning by specifying engine='python'.
  names = pd.read_csv('names.dmp',sep='\t\|\t', lineterminator='\t\|\n')
lca.py:126: ParserWarning: Falling back to the 'python' engine because the 'c' engine does not support regex separators (separators > 1 char and different from '\s+' are interpreted as regex); you can avoid this warning by specifying engine='python'.
  nodes = pd.read_csv('nodes.dmp',sep='\t\|\t', lineterminator='\t\|\n')
constructing tree....
start constructing tree....
processing rows 0 to 224839
processing rows 224839 to 449678
processing rows 449678 to 674517
processing rows 674517 to 899356
processing rows 899356 to 1124195
processing rows 1124195 to 1349034
processing rows 1349034 to 1573873
processing rows 1573873 to 1798712
processing rows 1798712 to 2023551
processing rows 2023551 to 2248390
processing rows 2248390 to 2473229
finish constructing tree....
There are multiple entries associated with the name. Do you mean Monera <bacteria> or Monera <archaea>?
Monera <archaea>
Monera <archaea>
19633    2157
Name: tax_id, dtype: int64
There are multiple entries associated with the name. Do you mean Bacteria <bacteria> or Bacteria <walking sticks>?
Bacteria <walking sticks
Invalid input. Do you mean Bacteria <bacteria> or Bacteria <walking sticks>?
Bacteria <walking sticks>
Bacteria <walking sticks>
936752    629395
Name: tax_id, dtype: int64
The lease common ancestor of Monera and Bacteria is biota or cellular organisms
```

### Visualization of the Phylogenetic classification

### References

https://www.ncbi.nlm.nih.gov/taxonomy 
https://evolution.berkeley.edu/evolibrary/article/phylogenetics_04
https://www.geeksforgeeks.org/lca-for-general-or-n-ary-trees-sparse-matrix-dp-approach-onlogn-ologn/

```

```
