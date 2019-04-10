#!/bin/bash
# Assembly pipeline for de novo assembly plus reference based scaffolding
# USAGE: ./assembly_pipeline.sh sampleID reads_folder 'species'
# Example: ./assembly_pipeline_SPAdes.sh 16T0014 reads 'Francisella tularensis' params.cfg 2>&1 |tee -a 16T0014.log
# the parameter for the run are read from a configuration file in the current directory that is named 'parameter.cfg' 


if [ $# -lt 3 ]; then
  echo "usage: $0 <sampleID reads_folder 'species'>"
  exit
fi


ID=$1
PROCESS=$2
Species_name=$3
configfile=`pwd`/parameter.cfg


if [[ ! -e $configfile ]]; 
then 
  echo "could not open config file $configfile, exit!"
fi


# get location of this script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


echo Load configuration file $configfile
source "$configfile"


FINISHED="FINISHED"
FQ=fastq
PROCDIR=$(pwd)

echo -n "START: "
date +"%d.%m.%Y %H:%M:%S"

if [[ -e ${FINISHED}/${ID} ]];
then
  echo "Run ${ID} already finished, exit."
  exit 1
fi

if [[ -z "$THREADS" ]] ;
then
    THREADS=1
fi


# check if "_R1_" or "_1" in paired-end
FILES1=$(ls ${PROCESS}/*${ID}*_R1_*.gz  2>/dev/null | sort | tr "\n" " ")
NF1=$(echo $FILES1 | awk '{print NF}')
VR1="_R1_"
VR2="_R2_"

if [[ $NF1 -lt 1 ]]
then
    FILES1=$(ls ${PROCESS}/*${ID}*_1.*.gz | sort | tr "\n" " ")
    NF1=$(echo $FILES1 | awk '{print NF}')
    if [ $NF1 -lt 1 ]
    then
        echo "file pattern must match *ID*_R1_*.fastq.gz or *ID*_1.fastq.gz"
        exit 1
    else
        VR1="_1."
        VR2="_2."
    fi
fi


file1=$(ls ${PROCESS}/*${ID}*${VR1}*.gz | head -1)
echo "first file:" $file1
ES=$(echo `basename $file1` | awk -F "_" '{print $1}')

if [[ "$ES" == "$ID" ]];
then
    echo "sampleID is " $ID
    FILES1=$(ls ${PROCESS}/${ID}_*${VR1}*.gz | sort | tr "\n" " ")
    FILES2=$(ls ${PROCESS}/${ID}_*${VR2}*.gz | sort | tr "\n" " ")
else
    ES=$(echo `basename $file1` | awk -F "_" '{print $2}')
    if [ "$ES" == "$ID" ]
    then
        echo "sampleID is " $ID
        FILES1=$(ls ${PROCESS}/*_${ID}_*${VR1}*.gz | sort | tr "\n" " ")
        FILES2=$(ls ${PROCESS}/*_${ID}_*${VR2}*.gz | sort | tr "\n" " ")
    else
        echo "could not find file with given sampleID $ID"
        exit 1
    fi
fi

base1=${ID}_1
base2=${ID}_2

NF1=$(echo $FILES1 | awk '{print NF}')
NF2=$(echo $FILES2 | awk '{print NF}')

if [[ $NF1 -lt 1 ]];
then
    echo "Could not find files, exit"
    exit 1
fi

echo "fwd reads: " $FILES1
echo "rev reads: " $FILES2
echo $NF1 $NF2


if [[ "$NF1" -ne "$NF2" ]];
then
  echo "ERROR: The number of files to process are differing! $NF1 -ne $NF2, please check the pipeline process!"
  echo "       Check the files with: ls ${PROCESS}/${ID}_*_R1_*.gz"
  echo "                             ls ${PROCESS}/${ID}_*_R2_*.gz"
  echo
  exit 1
fi


# download Genbank reference assemblies if they are not given in config file or if files do not exist
if [[ ! -z $REFERENCES ]] && [[ "$(ls -A $REFERENCES)" ]];
then
  echo "REFERENCES are in directory $REFERENCES"
else
    echo "Download genbank references"
    REFERENCES=${PROCDIR}"/Genbank_references"
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -P $REFERENCES
    echo "awk -v organism_name=\"$Species_name\" -v status=\"Complete Genome|Chromosome\" -v output=\"$REFERENCES\" -f $DIR/scripts/download_assemblies.awk $REFERENCES/assembly_summary.txt"
    awk -v organism_name="$Species_name" -v status="Complete Genome|Chromosome" -v output="$REFERENCES" -f $DIR/scripts/download_assemblies.awk $REFERENCES/assembly_summary.txt
    # save References directory to config file for next runs
    echo "REFERENCES=$REFERENCES" >>$configfile
fi


echo "------------------------------------"
date
echo "Uncompress and merge files from different runs"
if [[ ! -e $FQ ]]; then mkdir $FQ; fi
if [[ ! -e $FQ/${base1}.fastq ]]; 
then 
    echo  "gunzip -c $FILES1 >$FQ/${base1}.fastq";
    gunzip -c $FILES1 >$FQ/${base1}.fastq; 
fi
if [[ ! -e $FQ/${base2}.fastq ]]; 
then 
    echo "gunzip -c $FILES2 >$FQ/${base2}.fastq";
    gunzip -c $FILES2 >$FQ/${base2}.fastq; 
fi
echo "------------------------------------"
date
echo "Raw data quality check with FastQC"
if [[ ! -e fastqc ]]; then mkdir fastqc; fi

if [[ ! -e fastqc/${base1}_fastqc.zip ]]; 
then 
    fastqc -t $THREADS -f fastq -o fastqc $FQ/${base1}.fastq; 
fi

if [[ ! -e fastqc/${base2}_fastqc.zip ]]; 
    then 
    fastqc -t $THREADS -f fastq -o fastqc $FQ/${base2}.fastq; 
fi
echo "------------------------------------"

# remove fastqScreen 
#date
#echo "Raw data quality check with FastQ screen"
#if [ ! -e fastq_screen ]; then mkdir fastq_screen; fi
#if [ ! -e fastq_screen/${base1}_screen.txt ]; 
#then 
#    fastq_screen --conf fastq_screen.conf --threads $THREADS -o fastq_screen $FQ/${base1}.fastq; 
#fi

#if [ ! -e fastq_screen/${base2}_screen.txt ]; 
#    then 
#    fastq_screen --conf fastq_screen.conf --threads $THREADS -o fastq_screen $FQ/${base2}.fastq; 
#fi
#echo "------------------------------------"

date
echo "Kraken screen with mini-kraken DB"
if [[ ! -e mini_kraken ]]; then mkdir mini_kraken; fi

if [[ ! -e mini_kraken/${ES}.report.txt ]];
then
    kraken --fastq-input --db $minikrakenDB --paired --check-names --threads $THREADS $FQ/${base1}.fastq $FQ/${base2}.fastq >mini_kraken/${ES}.sequences.kraken
    kraken-translate --db $minikrakenDB mini_kraken/${ES}.sequences.kraken >mini_kraken/${ES}.sequences.labels
    kraken-report --db $minikrakenDB mini_kraken/${ES}.sequences.kraken >mini_kraken/${ES}.report.txt
fi
echo "------------------------------------"

date
echo "Adapter trimming & filtering with FlexBar 3.0.3"
if [[ ! -e flexbar ]]; then mkdir flexbar; fi

if [[ ! -e flexbar/${ES}.ec.fastq ]];
then
 flexbar --reads $FQ/${base1}.fastq --reads2 $FQ/${base2}.fastq --threads $THREADS --adapter-min-overlap 4 --adapters $DIR/configs/adapter.fasta --min-read-length 35 --adapter-trim-end ANY -f i1.8 -q 25 -i --target flexbar/${ES}.ec >flexbar/${ES}.ec.fastq
fi

echo "------------------------------------"

#date
#echo "Read error correction with Fiona"
# not working with large files (~11GB)
#if [ ! -e fiona ]; then mkdir fiona; fi

#if [ ! -e fiona/${ES}.ec.fastq ];
#then
#    fiona -g 1895994 -v -nt $THREADS flexbar/${ES}.ec.fastq fiona/${ES}.ec.fastq;
#fi
#echo "------------------------------------"


date
echo "Merge overlapping trimmed reads using Flash2"
if [[ ! -e flash2 ]]; then mkdir flash2; fi

if [[ ! -e flash2/${ES}.histogram ]];
then
  flash2 --interleaved -m 50 -r 300 -f 450 -s 100 -x 0.15 -t $THREADS --output-prefix flash2/${ES} flexbar/${ES}.ec.fastq;
fi
echo "------------------------------------"


#date
#echo "Compute de novo assembly using SPAdes --only-assembler --> runs only assembling (without read error correction)"
#if [ ! -e SPAdes ]; then mkdir SPAdes; fi
#if [ ! -e SPAdes/${ES} ]; then mkdir SPAdes/${ES}; fi

#if [ ! -e SPAdes/${ES}/contigs.fasta ];
#then
#  pe_reads=flash2_fiona/${ES}.notCombined.fastq
#  se_reads=flash2_fiona/${ES}.extendedFrags.fastq
#  spades.py --careful --only-assembler --12 ${pe_reads} -s ${se_reads} -o SPAdes/${ES} -t $THREADS
#  rm -R SPAdes/${ES}/split_input
#  rm -R SPAdes/${ES}/K*
#fi


date
echo "Compute de novo assembly using SPAdes --careful with error correction (default)"
if [[ ! -e SPAdes_err ]]; then mkdir SPAdes_err; fi
#if [ ! -e SPAdes_err/${ES} ]; then mkdir SPAdes_err/${ES}; fi

if [[ ! -e SPAdes_err/${ES}/contigs.fasta ]];
then
  pe_reads=flash2/${ES}.notCombined.fastq
  se_reads=flash2/${ES}.extendedFrags.fastq
  spades.py --careful --cov-cutoff auto --12 ${pe_reads} -s ${se_reads} -o SPAdes_err/${ES} -t $THREADS
  # --> error, run with flexbar output
  #spades.py --careful --12 flexbar/${ES}.ec.fastq -o SPAdes_err/${ES} -t $THREADS
  rm -R SPAdes_err/${ES}/split_input
  #rm -R SPAdes_err/${ES}/corrected
  #rm -R SPAdes_err/${ES}/K*
  #if [ -d SPAdes_err/${ES}/mismatch_corrector ];
  #then
    #rm -R SPAdes_err/${ES}/mismatch_corrector
  #fi
fi
echo "------------------------------------"

date
echo "Filter contigs by minimum coverage of $MIN_COV and length $MIN_LENGTH" 
if [[ ! -e ${PROCDIR}/SPAdes_err/${ES}/contigs_filtered.fasta ]];
then
    $DIR/scripts/filter_contigs.py -i ${PROCDIR}/SPAdes_err/${ES}/contigs.fasta -o ${PROCDIR}/SPAdes_err/${ES}/contigs_filtered.fasta -c $MIN_COV -l $MIN_LENGTH
fi

contig_assembly=${PROCDIR}/SPAdes_err/${ES}/contigs_filtered.fasta
echo "------------------------------------"

date
echo "Compute closest reference genome using MUMi value"
if [[ ! -e MUMi/${ES} ]]; then mkdir -p MUMi/${ES}; fi

if [[ -e MUMi/${ES}/${ES}_refgenome.txt ]];
then
  echo "read from MUMi/${ES}/${ES}_refgenome.txt"
  read -r refFasta_id < MUMi/${ES}/${ES}_refgenome.txt
fi

if [ ! $refFasta_id ];
then
    cd MUMi/${ES}
    refFasta_id=$(python $DIR/scripts/compute_mumi_multiproc.py -i $contig_assembly -g $REFERENCES -o ${ES}_mumi.txt -p ${THREADS})
    echo $refFasta_id >${ES}_refgenome.txt
    cd $PROCDIR
fi
echo "------------------------------------"


refFasta=$(ls ${REFERENCES}/*/${refFasta_id}.fna)
refGff=$(ls ${REFERENCES}/*/${refFasta_id}.gff)
refGbf=$(ls ${REFERENCES}/*/${refFasta_id}.gbff)

if [ ! -f $refFasta ];
then
    echo "Ref fasta is missing, exit!"
    exit 1
fi

if [ ! -f $refGff ];
then
    echo "Ref gff is missing, exit!"
    exit 1
fi

if [ ! -f $refGbf ];
then
    echo "Ref gbf is missing, exit!"
    exit 1
fi


date
echo "Map to closest reference genome $refFasta"

if [[ ! -e ref_mapping ]]; then mkdir ref_mapping; fi
mapping_bam=ref_mapping/${ES}.bam

if [[ ! -e qualimap/${ES}/qualimapReport.html ]];
then
    if [[ -e ${refFasta} ]]
    then
        if [[ ! -e ${refFasta}.1.bt2 ]];
        then
            ${bowtie2_dir}/bowtie2-build $refFasta $refFasta
        fi
        mapping_sam=ref_mapping/${ES}.sam
        mapping_bam=ref_mapping/${ES}.bam
        ${bowtie2_dir}/bowtie2 -x $refFasta --interleaved ${PROCDIR}/flexbar/${ES}.ec.fastq -S $mapping_sam -p $THREADS
        samtools sort $mapping_sam -o $mapping_bam
        rm $mapping_sam
    else
        echo "reference genome $refFasta does not exist, can not compute mapping"
    fi
fi
echo "------------------------------------"

date
echo "Mapping QC with qualimap"
if [[ ! -e qualimap ]]; then mkdir qualimap; fi
if [[ ! -e qualimap/${ES}/qualimapReport.html ]];
then
    if [ -e $mapping_bam ];
    then
        qualimap bamqc -bam $mapping_bam -outdir qualimap/${ES}
    fi
fi
echo "------------------------------------"


date
echo "Scaffolding with Reference genome"

echo "This program will check if PAGIT (Post Assembly Genome Improvement Toolkit) is installed correctly."
if [[ -z "$PAGIT_HOME" ]] ; then
   echo "Sorry, you have to source the installation file sourceme.pagit, see documentation.";
   exit ;
fi 

ABACAS_dir=${PROCDIR}/ABACAS/$ES
output_ABACAS=${ABACAS_dir}/scaffolds.fasta
mkdir -p $ABACAS_dir

if [[ ! -f $ABACAS_dir/$ES.fasta ]];
then
    echo "Running ABACAS:"
    cd $ABACAS_dir
    ln -sf $contig_assembly contigs.fasta
    ln -sf $refFasta
    perl $PAGIT_HOME/ABACAS/joinMultifasta.pl $refFasta Refsequence.union.fasta
    refFastaUnion=Refsequence.union.fasta
    perl $PAGIT_HOME/ABACAS/abacas.pl -t -r $refFastaUnion -q contigs.fasta -p nucmer -m -b -c -o $ES &> ../${ES}.out.abacas.txt
    perl $PAGIT_HOME/ABACAS/splitABACASunion.pl $refFasta $refFastaUnion ${ES}.fasta ${ES}.crunch ${ES}.tab
    cat Split.ABACAS.fasta ${ES}.contigsInbin.fas >scaffolds.fasta
fi

if [[ ! -f $output_ABACAS ]] ; then
   echo "Problem with ABACAS!";  
   exit 1
fi
cd $PROCDIR

echo "ABACAS ran successfully!"
echo "------------------------------------"



# IMAGE does not converge sometimes --> remove from pipeline !


echo "------------------------------------"
echo "Assembly improvement with Pilon"
echo "------------------------------------"

# Pilon automatically improves draft assemblies (and finds variation among strains, including large event detection)
# Pilon requires as input a FASTA file of the assembled genome along with a BAM file of aligned reads
# Pilon uses reads alignment analysis to identify inconsistencies between the input geonem and the evidence in the reads
# it then attempts to make inprovements to the input genome including
#   - single base sequences, small indels, larger idel or block substition events, gap filling and identification of local misassemblies
#     (optional opening of new gaps)
# outputs a FASTA file that contains an improved representation of the genome from the read data and an optional VCF file
# optionally produces tracks to be viewed in IGV or GenomeView and reports events suche as large collapsed repeat regions in its standard output_ABACAS


date
echo "Assembly improvement of ref-based scaffolds with Pilon"

Pilon_dir_scaffolds=${PROCDIR}/Pilon/${ES}
mkdir -p $Pilon_dir_scaffolds

scaffold_assembly=$output_ABACAS
assembly_sam=${Pilon_dir_scaffolds}/frags.sam
assembly_bam=${Pilon_dir_scaffolds}/frags.bam


if [[ ! -e ${Pilon_dir_scaffolds}/pilon.fasta ]] || [[ ! -s ${Pilon_dir_scaffolds}/pilon.fasta ]];
then
    if [[ ! -e ${Pilon_dir_scaffolds}/assembly.1.bt2 ]];
        then
        ${bowtie2_dir}/bowtie2-build $scaffold_assembly ${Pilon_dir_scaffolds}/assembly
    fi
    ${bowtie2_dir}/bowtie2 -x ${Pilon_dir_scaffolds}/assembly --interleaved ${PROCDIR}/flexbar/${ES}.ec.fastq -S $assembly_sam -p $THREADS
    samtools view -bS $assembly_sam | samtools sort -o ${Pilon_dir_scaffolds}/frags.bam
    samtools index ${Pilon_dir_scaffolds}/frags.bam
    java -Xmx16G -jar $PILON_JAR --genome $scaffold_assembly --frags $assembly_bam --outdir $Pilon_dir_scaffolds --changes --threads $THREADS
    rm ${assembly_sam}
    rm ${assembly_bam}
    rm ${Pilon_dir_scaffolds}/assembly*
fi

final_assembly=$Pilon_dir_scaffolds/pilon.fasta


echo "------------------------------------"
mapping_bam=${PROCDIR}/ref_mapping/${ES}.bam


date
echo "Assembly QC of contig assembly using QUAST"
echo "Reference-based assembly QC with QUAST and closest reference genome"
if [[ ! -e QUAST_contigs/${ES} ]]; then mkdir -p QUAST_contigs/${ES}; fi
if [[ ! -e QUAST_contigs/${ES}/report.tsv ]]; 
then
    quast.py -t $THREADS -R ${refFasta} -G ${refGff} --bam $mapping_bam -o QUAST_contigs/${ES} -L $contig_assembly
fi
cd $PROCDIR

echo "------------------------------------"

date
echo "Final Assembly QC with Quast"
echo "Reference-based assembly QC with QUAST and closest reference genome"

if [[ ! -e QUAST_final ]]; then mkdir QUAST_final; fi
if [[ ! -e QUAST_final/${ES}/report.tsv ]]; 
then
     # quast.py -t $THREADS -R ${refFasta} -G ${refGff} --bam $mapping_bam -o QUAST_final/${ES} -L --scaffolds $final_assembly
     # remove scaffolds option, because it is confusing in multiQC report
     quast.py -t $THREADS -R ${refFasta} -G ${refGff} --bam $mapping_bam -o QUAST_final/${ES} -L $final_assembly
fi
echo "------------------------------------"


date
echo "Run Prokka on contig assembly..."
Prokka_dir_contigs=${PROCDIR}/Prokka_contigs/${ES}
Prokka_result_contigs=${Prokka_dir_contigs}/${ES}.faa

if [[ ! -e ${Prokka_result_contigs} ]];
then
    if [[ ! -z "$Species_name" ]];
        then
            echo "VAR is not empty"
            stringarray=($Species_name)
            genus=${stringarray[0]}
            species=${stringarray[1]}
            # use genus for annotation if available
            prokka $contig_assembly --outdir $Prokka_dir_contigs --prefix $ES --genus "$genus" --species "$species" --strain $ES --kingdom 'Bacteria' --usegenus --force --cpus $THREADS
        else
            prokka $contig_assembly --outdir $Prokka_dir_contigs --prefix $ES --kingdom 'Bacteria' --strain $ES --force --cpus $THREADS
    fi
fi
echo "------------------------------------"

date 
echo "Run Prokka on final assembly..."

Prokka_dir_scaffolds=${PROCDIR}/Prokka_scaffolds/${ES}
Prokka_result_scaffolds=${Prokka_dir_scaffolds}/${ES}.faa

if [[ ! -e ${Prokka_result_scaffolds} ]];
then
    if [[ ! -z "$Species_name" ]]; 
        then
            echo "VAR is not empty"
            stringarray=($Species_name)
            genus=${stringarray[0]}
            species=${stringarray[1]}
            # use genus for annotation if available
            prokka $final_assembly --outdir $Prokka_dir_scaffolds --prefix $ES --genus "$genus" --species "$species" --strain $ES --kingdom 'Bacteria' --usegenus --force --cpus $THREADS
        else
            prokka $final_assembly --outdir $Prokka_dir_scaffolds --prefix $ES --kingdom 'Bacteria' --force --strain $ES --cpus $THREADS
    fi
fi
echo "------------------------------------"


mkdir -p $FINISHED

if [[ -e ${Prokka_result_scaffolds} ]];
then
  echo "Finished ${ES}." >$FINISHED/${ES}
fi 


echo "------------------------------------"
echo "------------------------------------"

echo -n "FINISHED. "
date +"%d.%m.%Y %H:%M:%S"
