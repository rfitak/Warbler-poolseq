# Enrichment analysis
Here we used the list of 93 genes ([outlier-genes.list](./data/outlier-genes.list)) to determine whether they were enriched for any specific gene ontology (GO) terms. Luckily, the entire great reed warbler proteome is available in [UniProt](https://www.uniprot.org) as proteome [UP000549775](https://www.uniprot.org/proteomes/UP000549775). We needed to first convert the NCBI gene IDs to UniProt IDs.

### Step 1: Convert NCBI to UniProt IDs
1. Download the table of the reference proteome in UniProt (https://www.uniprot.org/uniprotkb?query=proteome:UP000549775), ensure that the column "_Gene Names_" is selected.
   - this table contains the protein information for 13,628 reference proteins.
   - this table can be found here: [UP000549775.tsv](./data/UP000549775.tsv)
  2. Next, make a list of the UniProt IDs that match the 93 NCBI gene IDs
     - `cut -f1,5 UP000549775.tsv | grep -f outlier-genes.list | cut -f1 > outlier-UniProt.list`
     - this list can be found here: [outlier-UniProt.list](./data/outlier-UniProt.list)

### Step 2: Perform enrichment tests in PantherDB
Because the great reed warbler proteome in acessioned in UniProt, we are able to perform analyses using [PantherDB](https://pantherdb.org/). In PantherDB, the GO term enrichment test is referred to as a _Statistical overrepresentation test_.  There is no coding needed, so here are the steps used:
1. From the [PantherDB](https://pantherdb.org/) homepage, select the file [outlier-UniProt.list](./data/outlier-UniProt.list) or paste in the 93 UniProt IDs into the _**Enter IDs:**_ box.
2. For _**Select List Type:**_, select _ID's from Reference Proteome Genome_
3. Usign the dropdown menu, select _Acrocephalus_arundinaceus (ACRAR)_
4. Skip Section _**2. Select Organism**_
5. In Section _**3. Select Analysis**_, select _Statistical overrepresentation test_
6. In the resulting dropdown menu (_Choose annotation set_), select _PANTHER GO-Slim Molecular Function_ (or repeat for each different category).
7. Select _**Submit**_
8. In the next screen, select _Acrocephalus_arundinaceus_ACRAR_ as the reference proteome in the box **_Reference Proteome HMM List_**
9. A **_Selection Summary:**_ box should appear. Make sure the settings are correct (use the test settings below) then select **_Launch Aalysis_**
    - **_Test Type:_** Fisher's Exact
    - **_Correction:_** Calculate False Discovery Rate
  
The only signifcant results is displayed below:

Analysis Type:	PANTHER Overrepresentation Test (Released 20240807)
Annotation Version and Release Date:	PANTHER version 19.0 Released 2024-06-20
Analyzed List:	outlier-UniProt.list
Reference List:	Acrocephalus arundinaceus (ACRAR)
Test Type:	FISHER
Correction:	FDR
| PANTHER GO-Slim Molecular Function | Reference Genome | outlier-UniProt.list (93) | Expected | Fold Enrichment | P-value | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| triglyceride lipase activity (GO:0004806) | 13 | 4 | 0.09 | 45.09 | 1.39E-06 | 8.25E-04 |
| carboxylic ester hydrolase activity (GO:0052689) | 68 | 5 | 0.46 | 10.77 | 9.86E-05 | 2.93E-02 |
