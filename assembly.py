import random
import copy
import time
import sys
import os
import re
import operator
import math
import matplotlib.pyplot as plt
import networkx as nx

from utils import *
from dvoug import DVOUG

current_dir = os.path.dirname(__file__)
library_path1 = os.path.join(current_dir, "/home/lzz/Projects/DVOUG/vendor/kISS/build")  # 假设动态链接库在 build 目录下
library_path2 = os.path.join(current_dir, "/home/lzz/Projects/DVOUG/vendor/WFA2-lib/build")
sys.path.append(library_path1)
sys.path.append(library_path2)
import fm_index
from wfa2py import WFAlignerGapAffine, MemoryModel, AlignmentScope

# 设置初始值、步长和最大k-mer值
initial_kmer_size = 21
step_size = 20  # 每次增加的步长
max_kmer_size = 41

# 初始化所有序列数据
eDNAs_all = []
raw_eDNAs_all = []
eDNAs_all1 = read_fastq_sequences('se1_20X_2.fq')
eDNAs_all2 = read_fastq_sequences('se2_20X_2.fq')
raw_eDNAs_all.extend(eDNAs_all1)
raw_eDNAs_all.extend(eDNAs_all2)
eDNAs_all.extend(eDNAs_all1)
eDNAs_all.extend(eDNAs_all2)
min_k = 11
mid_k = 15

flag = False

# dbG = DVOUG()
# dbG.kmer_len = 221
# dbG.min_k = min_k
# dbG.mid_k = mid_k
# assembly_file = f'minia_assembly_k1.contigs.fa'
                
#                 # 将新的序列文件加载到eDNAs_all中
# eDNAs_all3 = fasta_to_dict(assembly_file)
# eDNAs_all.extend(list(eDNAs_all3.keys()))
# eDNAs_all.extend(list(eDNAs_all3.keys()))
# eDNAs_all.extend(list(eDNAs_all3.keys()))
# seq_to_fastq(eDNAs_all, 'seqfile.fq')
    
# # # # 构造bcalm命令
# command = (
#     f'./bcalm '
#     f'-in seqfile.fq '
#     f'-kmer-size 241 '
#     f'-abundance-min 1'
# )
    
# # # # 执行命令
# exit_status = os.system(command)
# with open('seqfile.unitigs.fa', 'r') as src, open("extend_seqfile.fa", 'w') as dest:
#     # 读取源文件内容并写入目标文件
#     dest.write(src.read())
# os.system(f'./minia -in extend_seqfile -kmer-size 241 -out minia_assembly_k241')
# dbG.reads = eDNAs_all
# dbG.wfa_aligner = WFAlignerGapAffine(1,1,1,AlignmentScope.Score,MemoryModel.MemoryMed)
# unitigs = fasta_to_dict("minia_assembly_k241.contigs.fa")
# dbG.unitigs_dict = unitigs
# dbG.unitigs = list(unitigs.keys())
# dbG.build_pseudogene()

# os.system('./kISS/build/kISS suffix_sort seqfile.fa -k 256 -t 4 --verbose')
# os.system('./kISS/build/kISS fmindex_build seqfile.fa -k 256 -t 4')
# dbG.fmi = fm_index.FMIndex_Uint32_KISS1(dbG.seq_end)
# dbG.fmi.load("/home/lzz/Projects/DVOUG/seqfile.fa.fmi")

# # dbG.ug_fmi.create_pseudogene(unitigs, work_dir+"ugfile.fa")
# os.system('./kISS/build/kISS suffix_sort ugfile.fa -k 256 -t 4 --verbose')
# os.system('./kISS/build/kISS fmindex_build ugfile.fa -k 256 -t 4')
# dbG.ug_fmi = fm_index.FMIndex_Uint32_KISS1(dbG.ug_end)
# dbG.ug_fmi.load("/home/lzz/Projects/DVOUG/ugfile.fa.fmi")

# dbG.unitigs_suffix_prefix()
        
# temp1 = time.time()
# dbG.extend_unitigs(1)
# temp2 = time.time()
# # dbG.merge_unitigs()

# G, contigs = all_contigs(dbG.new_unitigs, dbG.kmer_len)

# # 输出到指定文件
# extend_filename = "extend" + "_seqfile.fa"
# print_fasta_format_with_links_to_file(G, contigs,dbG.kmer_len, dbG.new_unitigs, extend_filename)
# os.system(f'./minia -in {extend_filename} -kmer-size 241 -out minia_assembly_k1')



# eDNAs_all1 = raw_eDNAs_all

# 循环从initial_kmer_size递增到max_kmer_size
for current_kmer_size in range(initial_kmer_size, max_kmer_size + 1, step_size):
    eDNAs_all1 = []
    eDNAs_all1.extend(raw_eDNAs_all)
    dbG = DVOUG()
    dbG.kmer_len = current_kmer_size
    if min_k <= 0 or flag:
        dbG.mid_threshold = 1
        # dbG.min_k = min_k
        dbG.min_k = current_kmer_size - 30
        if dbG.min_k < mid_k:
            dbG.min_k = mid_k
        dbG.mid_k = dbG.min_k
        # dbG.max_mis_num = math.ceil((current_kmer_size-min_k)/30)
        # dbG.min_k = mid_k
        # dbG.mid_k = mid_k
            # 计算新的文件名，`current_kmer_size - step_size`
        assembly_file = f'minia_assembly_k{current_kmer_size - step_size}.contigs.fa'
                
                # 将新的序列文件加载到eDNAs_all中
        eDNAs_all3 = fasta_to_dict(assembly_file)
        eDNAs_all1.extend(list(eDNAs_all3.keys()))
        eDNAs_all1.extend(list(eDNAs_all3.keys()))
        eDNAs_all1.extend(list(eDNAs_all3.keys()))
        # 将合并的序列写入新fastq文件
    else:
        dbG.min_k = min_k
        dbG.mid_k = mid_k
    seq_to_fastq(eDNAs_all1, 'seqfile.fq')
    command = (
        f'./vendor/Bloocoo '
        f'-file seqfile.fq '
        f'-kmer-size {current_kmer_size} '
        f'-abundance-min 2'
    )
    os.system(command)
    # 构造bcalm命令
    command = (
        f'./bcalm '
        f'-in seqfile_corrected.fastq '
        f'-kmer-size {current_kmer_size} '
        f'-abundance-min 2'
    )
    
    # 执行命令
    exit_status = os.system(command)
    if exit_status != 0:
        print(f"Error in bcalm execution for kmer-size {current_kmer_size}")
        break
    # with open('seqfile_corrected.unitigs.fa', 'r') as src, open('extend_seqfile.fa', 'w') as dest:
    #     # 读取源文件内容并写入目标文件
    #     dest.write(src.read())
    # os.system(f'./minia -in extend_seqfile.fa -kmer-size {current_kmer_size} -out minia_assembly_k{current_kmer_size}')

    # eDNAs_all_temp = fasta_to_dict(f'minia_assembly_k{current_kmer_size}.contigs.fa')
    # eDNAs_all1.extend(list(eDNAs_all_temp.keys()))
    # eDNAs_all1.extend(list(eDNAs_all_temp.keys()))
    # eDNAs_all1.extend(list(eDNAs_all_temp.keys()))
    dbG.reads = eDNAs_all1
    # seq_to_fastq(eDNAs_all1, 'seqfile.fq')
    # command = (
    #     f'./bcalm '
    #     f'-in seqfile.fq '
    #     f'-kmer-size {current_kmer_size} '
    #     f'-abundance-min 2'
    # )
    # os.system(command)
    dbG.wfa_aligner = WFAlignerGapAffine(1,1,1,AlignmentScope.Score,MemoryModel.MemoryMed)
    # unitigs = fasta_to_dict(f"minia_assembly_k{current_kmer_size}.contigs.fa")
    unitigs = fasta_to_dict("seqfile_corrected.unitigs.fa")
    dbG.unitigs_dict = unitigs
    dbG.unitigs = list(unitigs.keys())
    dbG.build_pseudogene()

    os.system('./vendor/kISS/build/kISS suffix_sort seqfile.fa -k 256 -t 4 --verbose')
    os.system('./vendor/kISS/build/kISS fmindex_build seqfile.fa -k 256 -t 4')
    dbG.fmi = fm_index.FMIndex_Uint32_KISS1(dbG.seq_end)
    dbG.fmi.load("/home/lzz/Projects/DVOUG/vendor/seqfile.fa.fmi")

    # dbG.ug_fmi.create_pseudogene(unitigs, work_dir+"ugfile.fa")
    os.system('./vendor/kISS/build/kISS suffix_sort ugfile.fa -k 256 -t 4 --verbose')
    os.system('./vendor/kISS/build/kISS fmindex_build ugfile.fa -k 256 -t 4')
    dbG.ug_fmi = fm_index.FMIndex_Uint32_KISS1(dbG.ug_end)
    dbG.ug_fmi.load("/home/lzz/Projects/DVOUG/vendor/ugfile.fa.fmi")

    dbG.unitigs_suffix_prefix()
        

    dbG.extend_unitigs()

    # dbG.merge_unitigs()

    G, contigs = all_contigs(dbG.new_unitigs, dbG.kmer_len)

    # 输出到指定文件
    extend_filename = "extend" + "_seqfile.fa"
    # extend_filename = "minia_assembly_k" + f'{current_kmer_size}' + ".contigs.fa"
    print_fasta_format_with_links_to_file(G, contigs,dbG.kmer_len, dbG.new_unitigs, extend_filename)
    # with open('seqfile.unitigs.fa', 'r') as src, open(extend_filename, 'w') as dest:
    #     # 读取源文件内容并写入目标文件
    #     dest.write(src.read())
    os.system(f'./vendor/minia -in {extend_filename} -kmer-size {current_kmer_size} -out minia_assembly_k{current_kmer_size}')
    flag = True