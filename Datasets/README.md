Here are the instructions to download and build up some datasets and reference databases.


## Simulated datasets

* sets of paired-end reads (setA2 and setB2_1)
* set of genomes 

### Sets of reads

The sets of simulated paired-end reads we used for experiments were designed ad hoc for a benchmarking analysis of metagenomic tools (Lindgreen, Adair & Gardner (2016)). They reproduce the size, complexity and characteristics of real metagenomic samples and moreover, in order to test the reliability of each classification, they include subsets of simulated random shuffled reads that should not be classified.

The datasets are available free for download at http://www.gardner-binflab.org/lindgreen-et-al-metagenomics-benchmark-data/). We selected the two sets of paired-end reads setA2 (setA2_1.fq, setA2_2.fq) and setB2 (setB2_1.fq, setB2_2.fq). Recall that due to size restrictions, each fastq file has been split in two parts (-0 and -1) that should be concatenated after download.

Before using them, we filtered the content of each fastq file by taking out the random shuffled reads and removing groups of reads associated with the phylum of Eukaryotes whose species provenance was not specified. In addition, we preferred to discard also a group of reads whose associated genomes are expiring entries in the NCBI database as being obsolete or incorrect. 
Thus, the set of reads considered for the positive control in setA2 (as well as setB2) has around 20 million of sequences, while the set of random shuffled reads considered for the negative control has around 5 million sequences. 

The list of reads we took into account for the positive control is partitioned in files title_setA2_part_.gz* and title_setB2_part_.gz*, respectively. To decompress them, one could use the commands

```sh
cat title_setA2_part_.gz* | zcat > title_setA2.txt
cat title_setB2_part_.gz* | zcat > title_setB2.txt
```

The list of random shuffled reads we took into account for the negative control is in title_setA2ran.txt.gz and title_setB2ran.txt.gz, respectively.

To build up the datasets starting from the list of reads (title_* .txt), one could use the script *filter_by_headers.pl* replacing the *fastq* and *headings* file names with their own file names. For instance, the set of reads considered for the positive control in setA2_1 can be obtained by running

```sh
perl filter_by_headers.pl > setA2_1.noRan.fq
```

### Reference database

The reference database comprises 930 genomes from 686 species belonging to 17 phyla as indicated in Lindgreen et al. https://www.nature.com/articles/srep19233 in the Supplementary Table S2. In Reference_database.csv we report the accession number and the lineage for any genome considered. 

To download genomes from their accession number, one could use the NCBI EFetch utility: try the following URL to download the nucleotide genome for CP000360.1
```sh
 wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000360.1&rettype=fasta&retmode=text"
 ```
 
## Real dataset

* set of paired-end reads (SRR1804065)
* set of genomes 


### Sets of reads

The set of real metagenome SRR1804065 is a tool sample of a female participant to the Human Microbiome Project (HMP), see https://www.ncbi.nlm.nih.gov/sra/SRR1804065. It is a paired-end read collection free to download using, for instance, the SRA toolkit utility *fastq-dump*.

Since the "ground truth" is not available for a real metagenome and a large number of reads is likely to belong to unknown species, we first filtered out the dataset by using BLAST.
The list of reads we took into account in SRR1804065 is partitioned in file title_SRR1804065.txt.gz*. 

To build up the dataset, one could decompress the list of examined reads (title_SRR1804065.txt)
```sh
cat title_SRR1804065.txt.gz* | zcat > title_SRR1804065.txt
 ```
and then filter out SRR1804065_1.fastq and SRR1804065_2.fastq using the script *filter_by_headers.pl* replacing the *fastq* and *headings* file names with their own file names.
```sh
perl filter_by_headers.pl > SRR1804065_1.filtr.fq
perl filter_by_headers.pl > SRR1804065_1.filtr.fq
 ```

### Reference database

The reference database we use for this experiment comprises 3,423 genomes from 1,499 species belonging to 42 phyla.
The list of their accession number is in Refs_SRR1804065.txt. As usual, to download genomes from their accession number, one could use the NCBI EFetch utility.
