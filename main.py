
import random
import copy
import time
import sys
import os
import re
import operator
import math
import argparse
from utils import *
from dvoug import DVOUG
current_dir = os.path.dirname(__file__)
library_path1 = os.path.join(current_dir, "/home/lzz/Projects/DVOUG/vendor/kISS/build")  # 假设动态链接库在 build 目录下
library_path2 = os.path.join(current_dir, "/home/lzz/Projects/DVOUG/vendor/WFA2-lib/build")
sys.path.append(library_path1)
sys.path.append(library_path2)
import fm_index
from wfa2py import WFAlignerGapAffine, MemoryModel, AlignmentScope

def main(eDNAs_all1, max_k, t2, t3, min_k, mid_k, t1, output_filename):
    dbG = DVOUG()
    dbG.kmer_len = max_k
    dbG.mid_threshold = t2
    dbG.min_threshold = t3
    dbG.min_k = min_k
    dbG.mid_k = mid_k

    seq_to_fastq(eDNAs_all1, 'seqfile.fq')

    command = (
        f'./bcalm '
        f'-in seqfile_corrected.fastq '
        f'-kmer-size {max_k} '
        f'-abundance-min {t1}'
    )
    
    exit_status = os.system(command)
    if exit_status != 0:
        print(f"Error in bcalm execution for kmer-size {max_k}")
        return

    dbG.wfa_aligner = WFAlignerGapAffine(1, 1, 1, AlignmentScope.Score, MemoryModel.MemoryMed)

    unitigs = fasta_to_dict("seqfile_corrected.unitigs.fa")
    dbG.unitigs_dict = unitigs
    dbG.unitigs = list(unitigs.keys())
    dbG.build_pseudogene()

    os.system('./vendor/kISS/build/kISS suffix_sort seqfile.fa -k 256 -t 4 --verbose')
    os.system('./vendor/kISS/build/kISS fmindex_build seqfile.fa -k 256 -t 4')
    dbG.fmi = fm_index.FMIndex_Uint32_KISS1(dbG.seq_end)
    dbG.fmi.load("/home/lzz/Projects/DVOUG/vendor/seqfile.fa.fmi")

    os.system('./vendor/kISS/build/kISS suffix_sort ugfile.fa -k 256 -t 4 --verbose')
    os.system('./vendor/kISS/build/kISS fmindex_build ugfile.fa -k 256 -t 4')
    dbG.ug_fmi = fm_index.FMIndex_Uint32_KISS1(dbG.ug_end)
    dbG.ug_fmi.load("/home/lzz/Projects/DVOUG/vendor/ugfile.fa.fmi")

    dbG.unitigs_suffix_prefix()
    dbG.extend_unitigs()

    G, contigs = all_contigs(dbG.new_unitigs, dbG.kmer_len)
    print_fasta_format_with_links_to_file(G, contigs, dbG.kmer_len, dbG.new_unitigs, output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DVOUG extension process.")
    parser.add_argument("--input", required=True, help="Input file for eDNAs_all1")
    parser.add_argument("--max_k", type=int, required=True, help="Max k value")
    parser.add_argument("--t2", type=int, required=True, help="Mid threshold value")
    parser.add_argument("--t3", type=int, required=True, help="Min threshold value")
    parser.add_argument("--min_k", type=int, required=True, help="Minimum k value")
    parser.add_argument("--mid_k", type=int, required=True, help="Mid k value")
    parser.add_argument("--t1", type=int, required=True, help="Minimum abundance for bcalm")
    parser.add_argument("--output", required=True, help="Output file name (FASTA)")

    args = parser.parse_args()

    main(
        eDNAs_all1=args.input,
        current_kmer_size=args.kmer_size,
        mid_threshold=args.mid_threshold,
        min_k=args.min_k,
        mid_k=args.mid_k,
        abundance_min=args.abundance_min,
        extend_filename=args.output
    )