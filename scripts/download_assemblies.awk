#!/usr/bin/awk -f
# awk -v organism_name="Francisella tularensis" -v status="Complete Genome|Chromosome" -v output="francisella_genomes" -f download_assemblies.awk assembly_summary.txt

BEGIN {FS = "\t"}
{

	if( tolower($8) ~ tolower(organism_name)".*" && tolower($12) ~ tolower(status)) {
	    print $1
	    output_dir=output"/"$1 
        fasta=$20 "/" $1 "_" $16 "_genomic.fna.gz"
	    gff=$20 "/" $1 "_" $16 "_genomic.gff.gz"
	    gbff=$20 "/" $1 "_" $16 "_genomic.gbff.gz"
	    system("mkdir -p " output_dir)
	    system("wget -nc " fasta " -P " output_dir " 2>/dev/null")
	    system("wget -nc " gff " -P " output_dir " 2>/dev/null")
	    system("wget -nc " gbff " -P " output_dir " 2>/dev/null")
	    system("gunzip -f " output_dir "/*.gz" " 2>/dev/null")
	}
}
END {
	print "organism_name: " organism_name
	print "assembly_level: " status
	print "output directory: " output
	#print ARGV[1]
	print "assembly summary file: " ARGV[1]
}
