#!/usr/bin/env python
#
# ==========================================================================================
# Kraken summary: Summarize abundance from Kraken report files for species and genus level
# ==========================================================================================

import sys
import argparse
import re
import math
import os
from collections import OrderedDict



def handle_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reports', '-r', type=argparse.FileType('r'), nargs='+', help="Kraken report files", required=True)
    parser.add_argument('--species_table', '-s', default="kraken_species_map_report_mqc.txt", help="Output summary table for species")
    parser.add_argument('--genus_table', '-g', default="kraken_genus_map_report_mqc.txt", help="Output summary table for genus")
    parser.add_argument('--species_min_sum', '-sm', type=int, default=20, help="Minimum sum of reads for a species in all samples (default: 20)")
    parser.add_argument('--genus_min_sum', '-gm', type=int, default=20, help="Minimum sum of reads for a genus in all samples (default: 20)") 
    parser.add_argument('--log10', action='store_true', help="Compute log10 of read counts")
    options = parser.parse_args()
    return options

def main(options):
    species_header = "# title: 'Kraken report (species level)'\n# description: 'Kraken report summary'\n# section: 'Custom Data File'\n# format: 'tsv'\n# plot_type: 'heatmap'\n# pconfig:\n#    id: 'samples'\n#    ylab: 'clades'"
    genus_header = "# title: 'Kraken report (genus level)'\n# description: 'kraken report summary'\n# section: 'Custom Data File'\n# format: 'tsv'\n# plot_type: 'heatmap'\n# pconfig:\n#    id: 'samples'\n#    ylab: 'clades'"  
    # hierarchy = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
    # hierarchy = ["S", "P", "C", "O", "F", "G", "S", "U"]
    species_candidates = re.compile(r"S|U")
    genus_candidates = re.compile(r"G|U")
    genus_reads = dict()
    species_reads = dict()
    species_list = []
    genus_list = []
    samples_list = []
    species_sum = dict()
    genus_sum = dict()

    for f in options.reports:
        sample = os.path.basename(f.name)
        sample = sample.replace(".report.txt","")
        samples_list.append(sample)
        for line in f:
            x = line.rstrip('\n')
            x_cells = x.split('\t')
            if re.search(species_candidates, x_cells[3]):
                species = x_cells[5]
                species = species.lstrip()
                species_reads[(sample, species)] = x_cells[1]
                if species not in species_list:
                    species_list.append(species)
                    species_sum[species] = int(x_cells[1])
                else:
                    species_sum[species] += int(x_cells[1])
            if re.search(genus_candidates, x_cells[3]):
                genus = x_cells[5]
                genus = genus.lstrip()
                genus_reads[(sample, genus)] = x_cells[1]
                if genus not in genus_list:
                    genus_list.append(genus)
                    genus_sum[genus] = int(x_cells[1])
                else:
                    genus_sum[genus] += int(x_cells[1])
    species_sum_sorted = OrderedDict(sorted(species_sum.items(), key=lambda t:t[1], reverse=True))
    species_list = species_sum_sorted.keys()
    genus_sum_sorted = OrderedDict(sorted(genus_sum.items(), key=lambda t:t[1], reverse=True))
    genus_list = genus_sum_sorted.keys()
    
    ###########################
    # write species map table #
    ###########################
    with open(options.species_table, 'w') as output_file:
        output_file.write('%s\n'%(species_header))
        output_file.write('%s'%('clade'))
        for s in samples_list:
            output_file.write('%s'%('\t' + s))
        for c in species_list:
            if species_sum[c] >= options.species_min_sum:
                output_file.write('\n')
                output_file.write('%s'%(c))
                for s in samples_list:
                    if species_reads.has_key((s,c)):
                        num_reads = species_reads[(s,c)]
                        if options.log10:
                            num_reads = math.log10( float(num_reads) )
                            output_file.write('\t%s'%(str(num_reads)) )
                    else:
                        output_file.write('\t%s'%('0'))
    output_file.close()
    #########################
    # write genus map table #
    #########################
    with open(options.genus_table, 'w') as output_file:
        output_file.write('%s\n'%(genus_header))
        output_file.write('%s'%('clade'))
        for s in samples_list:
            output_file.write('%s'%('\t' + s))
        for c in genus_list:
            if genus_sum[c] >= options.genus_min_sum:
                output_file.write('\n')
                output_file.write('%s'%(c))
                for s in samples_list:
                    if genus_reads.has_key((s,c)):
                        num_reads = genus_reads[(s,c)]
                        if options.log10:
                            num_reads = math.log10( float(num_reads) )
                        output_file.write('\t%s'%(str(num_reads)) )
                    else:
                        output_file.write('\t%s'%('0'))
    output_file.close()

if __name__ == '__main__':
    options = handle_args()
    main(options)
