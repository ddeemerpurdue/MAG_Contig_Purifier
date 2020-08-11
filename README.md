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
