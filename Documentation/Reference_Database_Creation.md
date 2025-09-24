# Reference Database Curation

These details were adapted from the Ecuador Bird Poop database created and shared by Dr. Sophie Picq and [Mike Allen](https://www.mikeallen-eco.com/news/2024/7/25/making-a-metabarcoding-reference-database-for-arthropods-part-1-downloading-sequences-1). Do not distribute outside of the lab without consulting Dr. Picq directly. Dr. Allen's blog is publicly available and I strongly recommend reading in full if this is new.

Create a reference database that has sequences from BOLD, NCBI, and MIDORI on metazoa and plants of the area (qza).

Use QIIME to then make a tree of sample ASVs once they are retrieved from DADA2.
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_single.qza \
  --output-dir mafft-fasttree-output

# turn into qzv to view
qiime tools export
  --input-path mafft-fastree-output/rooted_tree.qza \
  --output-path exported_tree

qiime tools export \
  --input-path mafft-fasttree-output/alignment.qza \
  --output-path exported_alignment

qiime tools export \
  --input-path mafft-fasttree-output/masked_alignment.qza \
  --output-path exported_masked_alignment
```
