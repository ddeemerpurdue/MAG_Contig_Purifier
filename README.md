## MAG_Contig_Purifier

### Purpose:

#### To take an existing set of bins from a metagenomic analysis and repatriate non-binning contigs, remove contigs that contaminate bins, and eventually merge/split bins using a combination of de-novo and reference-based methods.

## Outline
Starting point:
1. Raw assembly file(s)
2. Bin set
3. Bin identification file
- This is a tab-delimited file with two columns: bin_name and contig_name
- The bin_name should be in the format: "Bin.[number]"
- The contig_name should be in the format: ">[contig_name]"

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1a: Remove contamination from bins based on CATBAT annotation scores
# Step 1b: Add contigs into bins based on CATBAT annotation scores
# Step 2: Add contigs into bins based on ANI relations across samples
# Step 3: Merge bins together based on multiple lines of evidence:
a. Shared taxonomy
b. ANI patterns across samples
c. Non-overlapping gene content


