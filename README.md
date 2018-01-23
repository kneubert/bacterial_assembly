# bacterial_assembly
Assembly pipeline for baterial isolates using SPAdes

## Installation
Download the zip file or clone it with   
```sh
git clone https://github.com/kneubert/bacterial_assembly
```

All scripts need to be made executable:  
```sh
chmod a+x bacterial_assembly/\*.sh  
chmod a+x bacterial_assembly/scripts/\*  
```

Put the pipeline in your PATH, e.g. in your bashrc:  
export PATH=$PATH:SRC_PATH/bacterial_assembly-master 

## Example run
First check, if all the programs under **'Prerequisites'** are installed in your path.

### **1. Configuration**
As a first step create a folder e.g. 'my_project' somewhere in your working directory and create the configuration file **'parameter.cfg'**   
```sh
mkdir my_project; cd my_project  
```
The configuration file can look like this: [parameter.cfg](https://raw.githubusercontent.com/kneubert/bacterial_assembly/master/example/parameter.cfg):

\# [general parameter]   
\# number of threads that can be used by the pipeline      
**THREADS**=32   

\# [QC & preprocessing]   
\# the location of the minikraken DB needs to be given 
**minikrakenDB**=/group/ag_abi/kneubert/soft/Kraken/krakenDB/minikraken_20171019_8GB

\# [Assembly]   
\# directory that contains reference assemblies in subdirectories named my accession numbers:   
\# **assembly_accession** e.g. GCA_000008985.1   
\# ___ **fasta file** e.g. GCA_000008985.1_ASM898v1_genomic.fna   
\# ___ **gene annotation file** e.g. GCA_000008985.1_ASM898v1_genomic.gff   
\# ___ **genbank file** e.g. GCA_000008985.1_ASM898v1_genomic.gbff   
\# if 'REFERENCES' is not defined, Genbank references will be downloaded automatically for the given species   
**REFERENCES**=/group/ag_abi/kneubert/References 

\# PAGIT istallation directory   
**PAGIT_HOME**=/group/ag_abi/kneubert/soft/PAGIT

\# Pilon executable jar file   
**PILON_JAR**=/group/ag_abi/kneubert/soft/Pilon/pilon-1.22.jar   

\# the minimum coverage of filtered contigs   
**MIN_COV**=5   

\# the minimim length of filtered contigs   
**MIN_LENGTH**=500   

\# [Mapping]   
\# Bowtie2 directory   
**bowtie2_dir**=/group/ag_abi/kneubert/soft/bowtie2-2.3.3.1-linux-x86_64   

### **2. Run the assembly pipeline**
To run a single sample call the pipeline script with the sample-ID, read directory and species as parameter:
```sh
assembly_pipeline_SPAdes.sh  [sample-ID] [read direcory]  [species]   
assembly_pipeline_SPAdes.sh  16T0014 reads  'Francisella tularensis'   
```
It can be useful to write all outputs to a log file:   
```sh
assembly_pipeline_SPAdes.sh 16T0014 reads 'Francisella tularensis' 2>&1 |tee -a 16T0014-sub1M.log   
```
To run multiple samples, just create a bash script file like '**jobs**' and source it:   
```sh
assembly_pipeline_SPAdes.sh 16T0014-sub1M reads 'Francisella tularensis' 2>&1 |tee -a 16T0014-sub1M.log   
assembly_pipeline_SPAdes.sh 11T0315-sub1M reads 'Francisella tularensis' 2>&1 |tee -a 11T0315-sub1M.log   
assembly_pipeline_SPAdes.sh FSC237-sub1M reads 'Francisella tularensis' 2>&1 |tee -a FSC237-sub1M.log   

source jobs  
```
***mulitQC***  
After the runs have finished start the [multiQC script](https://raw.githubusercontent.com/kneubert/bacterial_assembly/master/multiqc.sh) in the project directory to summarize QC statistics before (preQC) and after the assembly (postQC).  
This script should produce three folders **preQC**, **postQC_contigs** and **postQC_scaffolds**, that contain the QC reports in HTML format for the raw data (FastQC, Kraken), the contig assembly and the scaffold assembly (QUAST, QualiMap, Prokka). The HTML-reports can be opened with any Browser that supports Javascript.

## Prerequisites
The following programs need to be installed:
### Quality check of raw data
* **FastQC** (https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
* **Kraken** with the mini-kraken DB (8 GB) (https://ccb.jhu.edu/software/kraken/)

### Adapter removal
* **Flexbar 3.0.3** (https://github.com/seqan/flexbar)

### Merge overlapping paired-end reads
* **FLASH2** (https://github.com/dstreett/FLASH2)

### De novo assembly (including read error correction)
* **SPAdes** with options --careful --cov-cutoff auto (http://cab.spbu.ru/software/spades/)

### Reference-based scaffolding of contigs
* **mummer** is needed to compute Maximum Unique Matches Index (MUMi) values between the de novo assembly and several Genbank reference assemblies (http://mummer.sourceforge.net/)
* **PAGIT**, the Post Assembly Genome Improvement Toolkit, is used for scaffolding (http://www.sanger.ac.uk/science/tools/pagit)
* **Pilon** is used for assembly improvement (https://github.com/broadinstitute/pilon)
* **Bowtie2** as a prerequisite for Pilon (http://bowtie-bio.sourceforge.net/bowtie2)
* **Samtools 1.3.1** is used to convert and index alignment files (https://sourceforge.net/projects/samtools/files/samtools/1.3.1/)

### Assembly QC
* **QUAST** to check the quality of the contig and the scaffold assembly separately (http://quast.sourceforge.net/quast.html)
* **QualiMap** to confirm mapping of reads to the reference assembly (http://qualimap.bioinfo.cipf.es/)
* **MultiQC** is used to combine quality statistics from different tools in a HTML report (http://multiqc.info/)

### Gene annotation
* **Prokka** is used to predict genes from the de novo contig assembly or scaffolds (https://github.com/tseemann/prokka)
