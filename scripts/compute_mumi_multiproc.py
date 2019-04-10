#!/usr/bin/python
# covnvert GenBank to FASTA nucleotides (*.gbff to *.ffn) 
# to obtain nucleotide sequences for all genes 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
import glob
import subprocess
import time

def handle_args():
    import argparse
    usage = ""
    usage += "Select the reference genome with the minimum MUMi value to the given draft assembly."

    parser = argparse.ArgumentParser(description = usage )
    #parser.print_help()
    parser.add_argument('-i', dest='assembly', help='assembly file [*.fna]')
    parser.add_argument('-g', dest='genomes', help='directory with reference genome assemblies [*.fna]')
    parser.add_argument('-o', dest='outfile', help='MUMi distances [*.tsv]')
    parser.add_argument('-p', dest='threads', help='Number of threads')
    args = parser.parse_args(sys.argv[1:])
  #  args = parser.parse_args()
    return args



def main(options):
    max_processes = options.threads
    genomes_dir = options.genomes
    output_filename = options.outfile
    info_filename = options.outfile
    assembly_filename = options.assembly
    output_dir = os.path.dirname(output_filename) 
    if os.path.isdir(output_dir):
        os.chdir(output_dir)
    #print "output dir: " +  output_dir

    script_path = os.path.dirname(os.path.realpath(__file__))

    ## process assembly file: concatenate sequences
    input_handle = open(assembly_filename, "r")
    
    output_assembly_tmp = output_dir + os.path.basename(assembly_filename) + ".tmp"
    #print "assembly tmp:" + output_assembly_tmp 
    output_handle = open(output_assembly_tmp, "w")
    records = list(SeqIO.parse(input_handle, "fasta"))
    seqlen_assembly = 0
    conseq = Seq("")
    nstring = 'N'*100
    n = len(records)-1
    start = 1
    end = 1
    
    for i in range(len(records)-1):
        conseq += records[i].seq
        conseq += Seq(nstring)
        seqlen_assembly += len(records[i].seq)
        seqlen_assembly += 100 # modified (length of inserted N's)
        end = start + len(records[i].seq) - 1
        start = end + 1
    
    conseq += records[n].seq
    seqlen_assembly += len(records[n].seq)
    desc = records[0].description + " (" + str(n+1) + " concatenated contigs)"
    
    processed_record = SeqRecord(conseq, id=records[0].id, description=desc)
    SeqIO.write(processed_record, output_handle, "fasta")

    end = start + len(records[n].seq) - 1
    #print records[n].id + "\t" + str(start) + "\t" + str(end)
    input_handle.close()
    output_handle.close()
    
    ## process reference genomes
    files = sorted(glob.glob(genomes_dir + "/**/*.fna"))
    n = 1
    min_mumi = 1
    closest_genome = ""
    MUM_num = 1

    processes = set()
    seq_length = {}
    proc_files = {}
    mum_files = {}
    
    # compute maximum unique matches using MUMer
    for file in files:
        #print "process " + file
        seqlen = 0
        conseq = Seq("")
        nstring = 'N'*100
        start = 1
        end = 1
        input_handle = open(file, "r")
        output_tmp = output_dir + os.path.basename(file) + ".tmp"
        output_handle = open(output_tmp, "w")
        records = list(SeqIO.parse(input_handle, "fasta"))
        seqlen = 0
        conseq = Seq("")
        nstring = 'N'*100
        n = len(records)-1
        start = 1
        end = 1
    
        for i in range(len(records)-1):
            conseq += records[i].seq
            conseq += Seq(nstring)
            seqlen += len(records[i].seq)
            seqlen += 100 # modified (length of inserted N's)
            end = start + len(records[i].seq) - 1
            start = end + 1
    
        conseq += records[n].seq
        seqlen += len(records[n].seq)
        desc = records[0].description + " (" + str(n+1) + " concatenated contigs)"
        processed_record = SeqRecord(conseq, id=records[0].id, description=desc)
        SeqIO.write(processed_record, output_handle, "fasta")
        end = start + len(records[n].seq) - 1
        input_handle.close()
        output_handle.close()
        MUM_file = "MUM_K12_CFT_l19_assembly_vs_" + str(MUM_num)
    
        mummer_command = "mummer -mum -b -c -l 19 " + output_assembly_tmp + " " + output_tmp + " >" + MUM_file
        command = "mummer"
        params = "-mum -b -c -l 19 " + output_assembly_tmp + " " + output_tmp
        seq_length[MUM_num] = seqlen
        proc_files[MUM_num] = file
        mum_files[MUM_num] = MUM_file
        MUM_num = MUM_num + 1
        
        log = open(MUM_file, 'a')
        #log.write('test')
        #log.flush()
     #   processes.add(subprocess.Popen([command, name]))
        processes.add(subprocess.Popen(['mummer -mum -b -c -l 19 %s %s' %(output_assembly_tmp, output_tmp)], stdout=log, shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])
        #print mummer_command
        #os.system(mummer_command)

    
        #mumi_command = "perl /group/ag_abi/kneubert/soft/MUMi/give_mumi.pl " + MUM_file + " -l1 " + str(seqlen_assembly) + " -l2 " + str(seqlen)
    # Check if all th child processes were closed
    
    for p in processes:
       if p.poll() is None:
           p.wait()
    
    #print "finished Mummer"
    
    min_mumi = 1
    info_handle = open(info_filename, "w")
    info_handle.write("file\tmumi\n")
    
    #print len(seq_length)
    
    # compute minimum MUMi value for all reference genomes
    for i in xrange(1,len(seq_length),1):
        mumi_command = "perl " + script_path + "/give_mumi.pl " + mum_files[i] + " -l1 " + str(seqlen_assembly) + " -l2 " + str(seq_length[i])
        var = os.popen(mumi_command).readlines()
        mumi = float(var[0])
        # print mumi
        
        #print "result: " + str(var[0])
        if mumi < min_mumi:
            min_mumi = mumi
            closest_genome = proc_files[i]
        info_handle.write(proc_files[i] + "\t" + str(mumi) + "\n")

    
    tmp_MUMfiles = sorted(glob.glob("MUM_*"))
    for out_tmp in tmp_MUMfiles:
         os.remove(out_tmp)
         
    tmp_files = sorted(glob.glob("*.tmp"))
    for out_tmp in tmp_files:
         os.remove(out_tmp)
         
    base = os.path.basename(closest_genome)
    closest_genome_name = os.path.splitext(base)[0]
    print closest_genome_name
            
    info_handle.close()
    # os.remove(output_assembly_tmp)

if __name__ == '__main__':
    options = handle_args()
    main(options)

