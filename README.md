# bacterial_assembly
Assembly pipeline for baterial isolates using SPAdes

## prerequisites
The following programs need to be installed
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
