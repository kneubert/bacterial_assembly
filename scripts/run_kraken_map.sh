#!/bin/bash
# USAGE:./run_kraken_map.sh kraken_results outputdir
input=$1
output=$2

source `pwd`/parameter.cfg

mkdir -p $output
file1=${output}/kraken_genus_map_report_mqc.txt
file2=${output}/kraken_species_map_report_mqc.txt

# get location of this script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# produce map report by kraken 
kraken-mpa-report --db $minikrakenDB $input/*sequences.kraken >$input/mpa_report.tsv

# filter entries for species level
grep "|s__" $input/mpa_report.tsv | awk -F "|" '{print $NF}' | sed 's/s__//g'  >$input/species_map_report.tsv

# filter entries for genus level
grep "|g__" $input/mpa_report.tsv | grep -v "s__" | awk -F "|" '{print $NF}' | sed 's/g__//g' >$input/genus_map_report.tsv

# filter entries on value and add sample IDs
filtered_genus_map=${input}/genus_map_report_filtered.tsv
filtered_species_map=${input}/species_map_report_filtered.tsv
samples=$(for file in  ${input}/*sequences.kraken; do echo `basename $file .sequences.kraken`; done)
echo "clade" $samples | tr " " "\t" >$filtered_genus_map
echo "clade" $samples | tr " " "\t" >$filtered_species_map

# filter only entries that have at least 20 in the sum of all samples 
awk '{ for(i=1; i<=NF;i++) j+=$i; if( j >= 20 ) {print $0;} j=0 }' $input/genus_map_report.tsv >>$filtered_genus_map

# filter only entries that have at least 5 in the sum of all samples 
awk '{ for(i=1; i<=NF;i++) j+=$i; if( j >= 20 ) {print $0;} j=0 }' $input/species_map_report.tsv >>$filtered_species_map


# Kraken genus report (heatmap)
#echo "# title: 'Kraken report (genus level)'" >$file1
#echo "# description: 'kraken-mpa-report --db $minikrakenDB mini_kraken/*sequences.kraken >mini_kraken/mpa_report.tsv'" >>$file1
#echo "# section: 'Custom Data File'" >>$file1
#echo "# format: 'tsv'" >>$file1
#echo "# plot_type: 'heatmap'" >>$file1
#echo "# pconfig:" >>$file1
#echo "#    id: 'samples'" >>$file1
#echo "#    ylab: 'clades'" >>$file1

#cat $filtered_genus_map >>$file1
R < $DIR/create_kraken_mqc_file.R --no-save

# Kraken species report (bargraph)
echo "# title: 'Kraken report (species level)'" >$file2
echo "# description: 'kraken-mpa-report --db $minikrakenDB mini_kraken/*sequences.kraken >mini_kraken/mpa_report.tsv'" >>$file2
echo "# section: 'Custom Data File'" >>$file2
echo "# format: 'tsv'" >>$file2
echo "# plot_type: 'bargraph'" >>$file2
echo "# pconfig:" >>$file2
echo "#    id: 'samples'" >>$file2
echo "#    ylab: 'clades'" >>$file2
cat $filtered_species_map >>$file2
