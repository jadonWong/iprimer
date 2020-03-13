#!/usr/bin/env python3
# -*-coding:utf-8-*-

"""
Usage:
    iprimer design [options] <genome.fasta> <loci.bed>
    iprimer (-h | --help)
    iprimer --version

Options:
-h,--help                        show this screen
--version                        show version
-f <int>,--flank <int>           size of flank sequence retained for loci [default: 10]
-r <int>,--primer_range <int>    range to find forward and reverse primers [default: 500]
-tm <float>,--min_TM <float>     the min TM of primer [default: 50.0]
-tM <float>,--max_TM <float>     thr max TM of primer [default: 68.0]
"""

from docopt import docopt
import primer3
import pandas as pd

if __name__ == '__main__':
    args = docopt(__doc__, version="2.0")
    #print(args)

    # make an index for genome
    faidx = {}
    count = 0
    seek_value = 0
    genome = args['<genome.fasta>']
    g = open(genome, 'r')
    while True:
        line = g.readline()
        if line:
            count += 1
            if count == 2:
                base_per_line = len(line) - 1
            seek_value += len(line)
            if line.startswith(">"):
                chr_name = line.split()[0].strip(">")
                faidx[chr_name] = seek_value
        else:
            break
    #print(faidx)
    g.close()
    
    # get the seq of each locus
    flank = int(args['--flank'])
    primer_range = int(args['--primer_range'])
    with open(args['<loci.bed>'], 'r') as bed:
        loci = bed.readlines()
    g = open(genome, 'r')
    for locus in loci:
        locus_info = locus.strip().split()
        chrome = locus_info[0]
        start = int(locus_info[1]) - flank - primer_range
        end = int(locus_info[2]) + flank + primer_range
        name = locus_info[3]
        product_min = int(locus_info[4])
        product_max = int(locus_info[5])

        #calculate the seek start and end, and read the sequence
        seek_start = faidx[chrome] + start // base_per_line * (base_per_line + 1) + start % base_per_line - 1
        seek_end = faidx[chrome] + end // base_per_line * (base_per_line + 1) + end % base_per_line
        g.seek(seek_start, 0)
        seq = "".join(g.read(seek_end - seek_start).split())
        #print(seq, len(seq))

        # design primer for each seq
        seq_args = {
            'SEQUENCE_ID':name,
            'SEQUENCE_TEMPLATE':seq,
            'SEQUENCE_INCLUDED_REGION':[0,len(seq)-1],
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':(1, primer_range, len(seq) - primer_range, primer_range,)
        }

        min_TM = float(args['--min_TM'])
        max_TM = float(args['--max_TM'])
        global_args = {
            'PRIMER_NUM_RETURN':10,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': min_TM,
            'PRIMER_MAX_TM': max_TM,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE':[[product_min, product_max]]
        }
    
        primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
        if primer3_result["PRIMER_PAIR_NUM_RETURNED"]:
            primer3_result_table_dict = {}
            index = []
            for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
                index.append("PRIMER_PAIR_" + str(i))
                primer_id = str(i)
                for key in primer3_result:
                    if primer_id in key:
                        info_tag = key.replace("_" + primer_id, "")
                        try:
                            primer3_result_table_dict[info_tag]
                        except:
                            primer3_result_table_dict[info_tag] = []
                        finally:
                            primer3_result_table_dict[info_tag].append(primer3_result[key])
        primer3_result_df = pd.DataFrame(primer3_result_table_dict, index=index).T
        primer3_result_df.to_csv(name + "_primer.csv", sep="\t")
    g.close()
