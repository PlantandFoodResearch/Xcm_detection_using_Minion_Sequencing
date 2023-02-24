Oxford Nanopore Technology (ONT) Sequencing for the detection of *Xanthomonas citri* subsp. *malvacearum* 
============================================================================================================

This analysis workflow was developed to process data generated using an ONT sequencing platform and analyse it for the presence of *Xanthomonas citri* subsp. *malvacearum* reads. 

**Outline**
============
1. Tool requirements and setting up your environment
2. Basecalling
3. Quality assessment
4. Concatenate reads into one file and filter
5. Metagenomic assembly 
6. Taxonomic assessment
   + Metamaps
   + Kraken
7. Race differentiation

Unless stated in a code chunk, run code from your chosen working directory.
Where the number of threads are specified these can modified as appropriate for your computing envronment.
  
**1. Tool requirements and setting up your environment**  
=======================================

This workflow requires a Linux environment and is dependent on the following tools:
+ [conda](https://docs.conda.io/en/latest/)
+ [docker](https://docs.docker.com/)
+ Guppy
+ Duplex tools
+ [NanoPlot](https://github.com/wdecoster/NanoPlot)
+ [NanoFilt](https://github.com/wdecoster/nanofilt)
+ [metaFlye](https://github.com/fenderglass/Flye)
+ [metamaps](https://github.com/DiltheyLab/MetaMaps) and/or
+ [kraken2](https://github.com/DerrickWood/kraken2/wiki)
+ [abricate](https://github.com/tseemann/abricate)
+ [krona](https://github.com/marbl/Krona/wiki)
+ [seqkit](https://bioinf.shenwei.me/seqkit/)

**Installation**

With the exception of guppy and duplex_tools the software dependencies may be managed through conda.

```{bash}
conda create -p Xcm_detection
conda activate Xcm_detection
conda install -c bioconda -c defaults nanoplot
conda install -c bioconda -c defaults nanofilt
conda install -c bioconda -c defaults flye
conda install -c bioconda -c defaults metamaps
conda install -c bioconda -c defaults kraken2
conda install -c bioconda -c defaults krona
conda install -c bioconda -c defaults seqkit
conda install -c bioconda -c defaults abricate
```

The version of Guppy suitable for your computing environment can be downloaded from the [Nanopore Community page](https://community.nanoporetech.com/downloads). Please follow the relevant installation instructions.

Duplex tools is available from the [Nanopore github page](https://github.com/nanoporetech/duplex-tools). 
Nanopore recommend running Duplex tools from an isolated virtual environment. Please refer to their [installation instructions](https://github.com/nanoporetech/duplex-tools#readme).


If you have trouble using metamaps with conda, metamaps can be set up on a [docker](https://docs.docker.com/) environment:

```{bash}
docker pull nanozoo/metamaps
```

Navigate to your working directory and create the necessary sub-directories
```{bash}
cd path/directory
mkdir fastq
mkdir fastq_duplex
mkdir nanoplot
mkdir flye
mkdir metamaps
mkdir kraken
mkdir abricate
```

**2. Basecalling**
========================================================================

Basecalling can be carried out either using the MinKnow software or with Guppy. 

Note here Qscore filtering during basecalling is disabled so quality assessment across the whole dataset can be done. Quality filtering is carried out at a subsequent step.
Use the relevant configuration file for your library preparation kit/flow cell combination to run Super Accurate basecalling. The protocol was developed using FLO-MIN106 flow cells/SQK-LSK110 library preparation kit and FLO-MIN114 flow cells/SQK-LSK114 library preparation kit.

### FLO-MIN106/SQK-LSK110

```{bash}
guppy_basecaller  --disable_qscore_filtering --input_path directory/reads.fast5 --save_path path/directory -c dna_r9.4.1_450bps_sup -x "cuda:0"
```

### FLO-MIN114/SQK-LSK114
Basecall simplex reads with read splitting enabled. 

```{bash}
guppy_basecaller  --disable_qscore_filtering --input_path directory/reads.fast5 --save_path path/directory -c dna_r10.4.1_e8.2_400bps_sup.cfg --do_read_splitting -x "cuda:0"
```
Identify a list of the template and complement reads using duplex_tools
```{bash}
duplex_tools pairs_from_summary sequencing_summary.txt path/directory
duplex_tools filter_pairs pair_ids.txt path/fastq_directory
```
Re-basecall using guppy_basecaller_duplex
```{bash}
guppy_basecaller_duplex -i <MinKNOW directory> -r -s duplex_calls -x 'cuda:0' -c dna_r10.4.1_e8.2_400bps_sup.cfg --chunks_per_runner 16 --duplex_pairing_mode from_pair_list --duplex_pairing_file pair_ids_filtered.txt
```

**3. Quality Assessment**
================================

Activate the conda environment to proceed with the pipeline.
```{bash}
conda activate Xcm_detection
```
The quality of the  run should be checked using an assessment tool such as Nanoplot (de Coster *et al.*, 2018). Nanoplot is available through the NanoPack package.

The location of the sequencing summary text file varies depending on whether the basecalling was post sequencing or during. You MUST modify this code chunk directing it to the sequencing summary file relevant to your dataset.

```{bash}
NanoPlot -t 10 --summary directory/fastq/sequencing_summary.txt -p SampleName -o nanoplot
```

If the run quality is low, preparing and sequencing another library is recommended.

**4. Concatenate all reads into one file and filter with NanoFilt**
==================================================================

Concatenate read files into a single fastq file and filter using a tool such as Nanofilt (de Coster *et al.*, 2018).

```{bash}
cat directory/*.fastq > directory/all_SampleName_raw.fastq
```
We recommend using different filtering criteria depending on the flowcells used.

### FLO-MIN106/SQK-LSK110
Remove read less than 1000bp or with a phred quality < 7 (default) 

```{bash}
cat directory/all_SampleName_raw.fastq | NanoFilt -l 1000 -s fastq/sequencing_summary.txt > directory/SampleName_filtered_l1000.fastq
```
### FLO-MIN114/SQK-LSK114
Remove read less than 1000bp or with a phred quality < 8 (default) 

```{bash}
cat directory/all_SampleName_raw.fastq | NanoFilt -l 1000 -q 8 -s fastq/sequencing_summary.txt > directory/SampleName_filtered_l1000.fastq
```

**5. Metagenomic Assembly**
================

Assemble the metagenome using a long-read metagenome assembler tool such as metaFlye (Kolmogorov *et al.*, 2020).

```{bash}
flye --nano-raw directory/SampleName_filtered_l1000.fastq --genome-size 6m --out-dir directory/flye --threads 20 --meta
```


**6. Taxonomic assessment**
==============================

### **MetaMap Analysis** 

The taxonomic profile of reads can be assigned using a tool such as metamaps (Dilthey *et al.*, 2019)

Here the miniSeq+H database (downloaded 16/10/2020) was used. 

```{bash}
metamaps  mapDirectly -t 40 --all -r directory/metamaps/databases/databases/miniSeq+H/DB.fa -q directory/SampleName_filtered_l1000.fastq -o directory/metamaps/classification_results
```

```{bash}
metamaps classify --mappings directory/metamaps/classification_results --DB directory/metamaps/databases/databases/miniSeq+H -t 40
```
Filter taxonomic assignments for a minion of 80% percent identity
```{bash}
perl path/MetaMaps-master/util/filterLowIdentityEntities.pl --DB miniSeq+H --mappings directory/metamaps/classification_results --identityThreshold 0.8
Rscript plotMappingSummary.R directory/metamaps/classification_results_filt80
```

Visualise using krona (Ondov *et al.*, 2011).

```{bash}
ktImportTaxonomy -i -m 2 -o directory/metamaps/SampleName.html directory/metamaps/classification_results_filt80.EM.reads2Taxon.krona
```

Extract reads classified as Xcm for further analysis using seqkit (Shen *et al.*, 2016).
```{bash}
awk '{if ($2 == "86040") print $0;}' metamaps/classification_results.EM.reads2Taxon | cut -f 1 > Xcm1.ids
wc -l Xcm1.ids
```

```{bash}
awk '{if ($2 == "1118965") print $0;}' metamaps/classification_results.EM.reads2Taxon | cut -f 1 > Xcm2.ids
wc -l Xcm2.ids
``` 

```{bash}
awk '{if ($2 == "1127439") print $0;}' metamaps/classification_results.EM.reads2Taxon | cut -f 1 > Xcm3.ids
wc -l Xcm3.ids
```

```{bash}
awk '{if ($2 == "1220027") print $0;}' metamaps/classification_results.EM.reads2Taxon | cut -f 1 > Xcm4.ids
wc -l Xcm4.ids
```

```{bash}
awk '{if ($2 == "1220028") print $0;}' metamaps/classification_results.EM.reads2Taxon | cut -f 1 > Xcm5.ids
wc -l Xcm5.ids
```

Combine the individual files into one file containing the read names for Xcm reads.
```{bash}
cat Xcm1.ids Xcm2.ids Xcm3.ids Xcm4.ids Xcm5.ids > Xcm_all.ids
```

Double check that the total number of reads in Xcm_all.ids matches the sum of the individual files from above:
```{bash}
wc -l Xcm_all.ids
```
Extract the Xcm reads from filtered data using seqkit
```{bash}
seqkit grep -f Xcm_all.ids directory/SampleName_filtered_l1000.fastq -o SampleName_Xcm_reads_metamaps.fastq
```

### **Kraken**

Alternatively Kraken2 (Wood *et al.*, 2019) can be used to analyse the raw reads or the metaFlye assembly.
Details on building a kraken2 database from genomes in RefSeq can be found at:
https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

Kraken2 can be run on either the metaflye assembly or the filtered reads.

Kraken2 on the metaflye assembly:
```{bash}
kraken2 --db directory/kraken2db/Refseq91/ --threads 20 --unclassified-out directory/kraken/SampleName_kraken2fl.unclassified.fasta --report directory/kraken/SampleName_kraken2fl.report --output directory/kraken/SampleName_kraken2fl.txt directory/flye/assembly.fasta
```
Kraken2 on the filtered reads:
```{bash}
kraken2 --db directory/kraken2db/Refseq91/ --threads 20 --unclassified-out directory/kraken/SampleName_kraken2fl.unclassified.fasta --report directory/kraken/SampleName_kraken2fl.report --output directory/kraken/SampleName_kraken2fl.txt directory/SampleName_filtered_l1000.fastq
```

Visualise using krona

```{bash}
ktImportTaxonomy -m 3 -t 5 -o directory/kraken/SampleName_kraken2f1.html directory/kraken/SampleName_kraken2fl.report
```

**7. Race differentiation**
================================

Search for race specific genes using Abricate (Seemann). A database has to be set up for this to work.

Make a database using the customized fasta file Race18_abricate.fasta provided here.


```{bash}
cd directory/abricate/db
cp directory/Race18_abricate.fasta sequences
makeblastdb -in sequences -title race_specific -dbtype nucl -hash_index
```

Use abricate to search for the sequences in Race18_abricate.fasta

```{bash}
abricate --db race_specific flye/assembly.fasta
```

**Acknowledgements**
================
This project is supported by the Grains Research and Development Corporation, through funding from the Australian Government Department of Agriculture, Fisheries & Forestry, as part of its Rural R&D for Profit program and along with Cotton Research and Development Corporation, Hort Innovation Australia, Wine Australia, Sugar Research Australia and Forest and Wood Products Australia.

References
=============

1.	De Coster W, D’Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-9.
2.	Dilthey AT, Jain C, Koren S, Phillippy AM. Strain-level metagenomic assignment and compositional estimation for long reads with MetaMaps. Nature Communications. 2019;10(1):3066.
3.	Kolmogorov M, Bickhart DM, Behsaz B, Gurevich A, Rayko M, Shin SB, et al. metaFlye: scalable long-read metagenome assembly using repeat graphs. Nature Methods. 2020;17(11):1103-10.
4.	Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011;12(1):385.
5.	Seemann T, Abricate, Github https://github.com/tseemann/abricate
6.	Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PloS one. 2016;11(10):e0163962.
7.	Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biology. 2019;20(1):257.
