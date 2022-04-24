# cgMLST-Klebsiella-oxytoca

cgMLST to Klebsiella oxytoca

The purpose of this repository is to describe how we created cgMLST for Klebsiella oxytoca. This scheme was created using the ChewBBACA pipeline (link below). 

chewBBACA
To download ChewBBACA access the link:https://github.com/B-UMMI/chewBBACA_tutorial

Workflow used to create the scheme
* Step 1: Scheme Creation of cgMLST
* Step 2: Allele calling
* Step 3: Scheme Validation (Allele call) of cgMLST
* Step 4: Extracting the Matrix target genes
* Step 5: Minimum Spanning Tree (MST)
* Step 6: Graphical evaluation of the scheme


Softwares and Downloads (Main dependencies)
* BLAST 2.12.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* Prodigal 2.6.3 https://github.com/hyattpd/prodigal/releases/ or above

# Step 1: Schema Creation
Selection of complete genomes for schema creation
As of March 3, 2022, 22 K. oxytoca complete genome sequences were available at the NCBI (National Center for Biotechnology Information) genome sequence repository (https://www.ncbi.nlm.nih.gov/assembly).
A list of all complete genomes used to create this protocol can be found in the file Complete_Genomes.xlsx.
Among the 22 complete genomes, K. oxytoca FDAARGOS_500 reference genome (GCF_003812925.1) was used by Prodigal algorithm as reference to recognize coding sequences (CDs). Prodigal generated the FDAARGOS_500.trn file at this step.
The FDAARGOS_500 genome was then removed from further analysis.
# Command:
```
#create schema
chewBBACA.py CreateSchema -i Complete_Genomes --cpu 46 -o schema_seed --ptf FDAARGOS_500.trn
```
The above command uses 46 CPU and creates a preliminary scheme (wgMLST) in the schema_seed folder using the trained product folder FDAARGOS_500.trn that was generated using the reference genome FDAARGOS_500 (GCF_003812925.1) and 21 complete genome sequences. The wgMLST scheme generated contained 9735 loci based on the 21 complete genomes.

# Step 2: Allele calling
In this step the allele calling is performed using the resulting set of loci determined in step 1.
# Command:
```
#run allelecall
chewBBACA.py AlleleCall -i Complete_Genomes -g schema_seed/ -o results_cg --cpu 46 --ptf FDAARGOS_500.trn
```
The allele calling used the default BLAST Score Ratio (BSR) threshold of 0.6.

# Step 2.1: Paralog detection
In this step genes considered paralogous from result of the allelecall (see above) are removed
# Command:
```
#run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_20220319T103614/results_alleles.tsv -g results_cg/results_20220319T103614/RepeatedLoci.txt -o alleleCallMatrix_cg
```
In this step, 67 genes were identified as possible paralogs and were removed from further analysis. 

# Step 2.2: Genome Quality Control
In this step we define a Threshold for the scheme that limits the loss of loci targets defined in the previous steps per genome and excludes genomes considered to be of low quality due to significant loci absence.
With this analysis we define the percentage of loci that will constitute the scheme based on how many targets we want to keep in this phase. For example, 100%, 99.5%, 99% and 95% of the loci may present in the set of high-quality genomes. This is one of the main steps in defining the cgMLST scheme targets.
# Command:
```
#run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```
In the Threshold 65 a set of 4085 loci were found to be present in all the analyzed complete genomes. 
In this stage we chose the loci present in 100% of the complete genomes and the Threshold 65 that limited the loss of the loci in the genomes. In this Threshold (65) 1 complete genome were removed due to loss of loci.
# Command:
```
#run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_65 --g GenomeRemoved65thr.txt  --t 1
```
This script selects all * loci * present in the selected * Threshold *. The value * t * is associated with the percentage of * loci * we want to be present in the set of genomes, for example: * t 1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_65 folder receives the result of the allelic profile for each of the 4085 candidate * loci * (allelic profile matrix). The file in this folder (cgMLST.tsv) contains the allelic profile of 4085 selected * loci * and will be used to create the core gene list. 

# Step 2.3: Creating the core gene list
This command selects all target genes from the "cgMLST.tsv" spreadsheet.
```
#list
head -n 1 cgMLST.tsv > Genes_100%_Core_65.txt
```
This step generated the file Genes_100%_Core_65.txt. This list needs to be transposed so that each core gene name is reported in a single line:
# Command:
```
#transpose table
datamash -W transpose < Genes_100%_Core_65.txt > Genes_Core_Al.txt 
```
This step generated the file > Genes_Core_Al.txt
You can see the list file with 4085 target genes at Genes_Core_Al.txt and for the subsequent steps we added the full path to each locus fasta file.
This list Genes_Core_Al.txt was then modified so that each name was preceeded by schema_seed:
# Command:
```
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("schema_seed/"$1)}' > list_genes_core.txt
```
# Step 3: Scheme Validation (Allele calling)
For the validation step we selected 386 drafts K. oxytoca genomes that were publicly available in RefSeq (https://www.ncbi.nlm.nih.gov/assembly) in March 2022. The list of all the draft genomes used can be found in the file Genomes_Validation.xlsx.
We this set of genomes (386 drafts genomes) we repeated the allele call step using only the 4085 candidate target genes.
# Command:
```
chewBBACA.py AlleleCall -i ../unfinished-genome/ --gl list_genes_core.txt -o ../results_all --cpu 46  -g ../schema_seed/schema_seed --ptf FDAARGOS_500.trn
```
The folder unfinished-genome contains the 386 validation drafts genomes used for validation of the scheme.

# Step 3.1: Concatenate the allelic profiles
The purpose of this step is to concatenate the matrix of the loci that defined the scheme and the matrix of the loci from the validation genomes. Thus, to concatenate the allelic profile matrix obtained from the creation of the scheme cgMLST_65/cgMLST.tsv with the matrix obtained for the validation genomes results_all/ results_20220319T114751/results_alleles.tsv. The following command was used:
# Command:
```
#create header
head -n 1 cgMLST_65/cgMLST.tsv > cgMLST_all.tsv
```
cgMLST_all.tsv file (cgMLST_all.tsv) contains the allelic profile of the 386 drafts genomes.

# Step 3.2: Evaluation of genome quality
After concatenation, we used the TestGenomeQuality to assess the impact of each validation genome on candidate loci in order to exclude low quality validation genomes. In this step you may need to define a new Threshold, as well as a new value of the parameter t, because loci that remain after the filters are the ones that will constituted the final scheme.
# Command:
```
chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```
# Step 4: Extracting the Matrix loci
At this step we chose loci present in 99% (--t 0.99) of the validation genomes and the Threshold 145 to limit the loss of the loci in the genomes.
In Threshold 145 a set of 3356 loci were found to be present in 99% the validation genomes.
We created another (removedGenomes145thr.txt) file with only the genomes removed at Threshold (145). 
Using Threshold (145) only 6 draft genomes were removed due to absence of loci targets.
 
# Command:
```
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_145 --t 0.99 --g removedGenomes145thr.txt 
```
This script selects loci and genomes that remained in the Threshold 145 and excludes the validation genomes and loci that were excluded with this Threshold.
The folder with the output file can be found at: cgMLST_145. This folder contains four files "cgMLST.tsv; cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv".
The cgMLST targets can be found at: cgMLST_145/cgMLSTschema.txt It contains the list of 3356 genes in the core genome defined as targets for this scheme. 

# Step 5: Minimum Spanning Tree
For the visualization of results, minimum spanning trees were buitl. Based on the allelic profiles obtained by the cgMLST scheme for each of the 400 genomes minimum spanning trees (MST) were constructed using the software GrapeTree (version 2.1.0) (https://github.com/achtman-lab/GrapeTree/releases) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection. The cgMLST_145/cgMLST.tsv file contains the allelic profile of the 400 genomes typed by cgMLST.

# Step 6: Graphical evaluation of the scheme
To assess the variability of the gene targets of cgMLST as well explore and evaluate the type and extent of allelic variation detected in each of the chosen loci. We run this script and graphically visualize the data in a series of html files.
# Command:
```
chewBBACA.py SchemaEvaluator -i schema_seed/ -o rms --cpu 6
```
