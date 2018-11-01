#!/bin/bash
# compute QC before (preQC) and after (postQC) assembly using multiQC
# execute inside working directory

# get location of this script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


# preQC (FastQC, minikraken)
$DIR/scripts/kraken_summary.py --reports mini_kraken/*report.txt --log10 --species_table mini_kraken/kraken_species_map_report_mqc.txt --genus_table mini_kraken/kraken_genus_map_report_mqc.txt
multiqc -n preQC -f -o preQC -c $DIR/configs/multiqc.config fastqc mini_kraken

# postQC (QUAST, qualimap, Prokka)
multiqc -n postQC_contigs -f -o postQC_contigs -c $DIR/configs/multiqc.config QUAST_contigs qualimap Prokka_contigs
multiqc -n postQC_scaffolds -f -o postQC_scaffolds -c $DIR/configs/multiqc.config QUAST_final qualimap Prokka_scaffolds
