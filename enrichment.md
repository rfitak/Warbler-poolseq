# Enrichment analysis
Here we used the list of 93 genes ([outlier-genes.list](./data/outlier-genes.list)) to determine whether they were enriched for any specific gene ontology (GO) terms. Luckily, the entire great reed warbler proteome is available in [UniProt](https://www.uniprot.org) as proteome [UP000549775](https://www.uniprot.org/proteomes/UP000549775). We needed to first convert the NCBI gene IDs to UniProt IDs:
1. Download the table of the reference proteome in UniProt (https://www.uniprot.org/uniprotkb?query=proteome:UP000549775), ensure that the column "_Gene Names_" is selected.
   - this table contains the protein information for 13,628 reference proteins.
   - this table can be found here: [UP000549775.tsv](./data/UP000549775.tsv)
  2. Next, make a list of the UniProt IDs that match the 93 NCBI gene IDs
     - `cut -f1,5 UP000549775.tsv | grep -f outlier-genes.list | cut -f1 > outlier-UniProt.list`
     - this list can be found here: [outlier-UniProt.list](./data/outlier-UniProt.list)
