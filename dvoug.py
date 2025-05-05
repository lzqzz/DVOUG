import re
import time
import sys
import os
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from Bio import Seq, SeqIO, SeqRecord
from utils import *
import networkx as nx
from itertools import islice
from collections import Counter
from bisect import bisect_right
from Unitigs_Info import UnitigInfo

current_dir = os.path.dirname(__file__)
library_path = os.path.join(current_dir, "/home/lzz/Projects/DVOUG/vendor/kISS/build")  # 假设动态链接库在 build 目录下
sys.path.append(library_path)

import fm_index

class DVOUG:
    def __init__(self, data=bytes(), kmer_len=21):
        self.data = data
        self.kmers = {}
        self.kmer_len = kmer_len
        self.min_k = 11
        self.mid_k = 15
        # self.mid_k1 = 121
        # self.mid_threshold1 = 2
        self.check_k = 2
        self.min_threshold = 2
        self.mid_threshold = 1
        self.max_gap_num = 5
        self.max_mis_num = 1
        self.gap_len = 100
        self.reads = []
        self.sequence_info = []
        self.seq_start = set()
        self.seq_end = set()
        self.ug_start = set()
        self.ug_end = set()
        self.unitigs_info1 = []
        self.unitigs_info = {}
        self.position_to_read = []
        self.position_to_ug = []
        self.break_num = 0
        self.unitigs = []
        self.unitigs_dict = {}
        self.new_unitigs = {}
        self.unitigs_ex_info = {}
        self.unitigs_prefix = set()
        self.unitigs_suffix = set()
        self.rc_unitigs_prefix = []
        self.rc_unitigs_suffix = []
        self.pseudogene = ''
        self.pathsA = []
        self.pathAB = []
        self.fmi = None
        self.ug_fmi = None
        self.wfa_aligner = None
        self.t1 = 0
        self.t2 = 0
        self.t3 = 0
        self.ug_break_point = {}
        self.ex_ug_break_point = {}


        self.set_max_kmer_num = False
        self.max_path_num = 1000
        

    def open_fastq(self, file):
        f = open(file, "r")
        matchLineB = re.compile('^(\+)')
        last_line = f.readline()
        line = f.readline()
        reads = 0
        while line.strip():
            if matchLineB.match(line):
                self.add_seq(last_line.strip())
                reads = reads + 1
                if reads > 1000:
                    #print('.', end='')
                    reads = 0
            else:
                last_line = line
            line = f.readline()

    def open_fastq_num(self, file, num =12000 * 100):
        f = open(file, "r")
        matchLineB = re.compile('^(\+)')
        last_line = f.readline()
        line = f.readline()
        reads = 0
        read_num = 1
        while line.strip() and read_num <= num:
            if matchLineB.match(line):
                self.add_seq(last_line.strip())
                reads = reads + 1
                if reads > 1000:
#                    print('.', end='')
                    reads = 0
            else:
                last_line = line
            line = f.readline()
            read_num = read_num + 1

    def open_dump(self, file):
        f = open(file, "r")
        line = f.readline()
        kmer_num = 0

        #Check K-mer length
        while (line.strip()):
            arr = line.split('\t')
            # print(arr)
            self.add_kmer(arr[0], int(arr[1]))
            rev_seq = DNA_rev_complement(arr[0])
            if rev_seq != arr[0]:
                self.add_kmer(rev_seq, int(arr[1]))
            kmer_num = kmer_num + 1
            if kmer_num > 10000000:
#                print('.', end='')
                kmer_num = 0
            line = f.readline()

    def open_fasta(self, file):
        print('Opening fasta file')
        f = open(file, "r")
        matchLineA = re.compile('^(>)')
        line = f.readline()
        seq_num = 0
        while (line.strip()):
            if not matchLineA.match(line):
                self.add_seq(line.strip())
            seq_num = seq_num + 1
            if seq_num > 1000:
#                print('.', end='')
                seq_num = 0
            line = f.readline()

    def open_row_seq_file(self, file):
        f = open(file, "r")
        line = f.readline()
        while (line.strip()):
            self.add_seq(line.strip())
            line = f.readline()

    def add_seqs(self, seqs):
        for seq in seqs:
            self.add_seq(seq)

    def add_seq(self, str):
        if len(str) >= self.kmer_len:
            i = 0
            kmstr = ''
            while i <= len(str) - self.kmer_len:
                kmstr = str[i:i + self.kmer_len]
                if not ("N" in kmstr):
                    self.add_kmer(kmstr)
                    self.add_kmer(DNA_rev_complement(kmstr))
                else:
                    if kmstr.count("N") == 1:
                        self.add_kmer(kmstr.replace("N", "A"))
                        self.add_kmer(kmstr.replace("N", "T"))
                        self.add_kmer(kmstr.replace("N", "G"))
                        self.add_kmer(kmstr.replace("N", "C"))
                        self.add_kmer(Seq.reverse_complement(kmstr.replace("N", "A")))
                        self.add_kmer(Seq.reverse_complement(kmstr.replace("N", "T")))
                        self.add_kmer(Seq.reverse_complement(kmstr.replace("N", "G")))
                        self.add_kmer(Seq.reverse_complement(kmstr.replace("N", "C")))
                i = i + 1
    
    def add_kmer(self, kmer, num=1):
        # comp_kmer = self.compressDNA(kmer)
        if kmer in self.kmers:
            self.kmers[kmer] = self.kmers[kmer] + num
        else:
            self.kmers[kmer] = num
    
    def build_pseudogene(self):
        pseudogene = ''
        unitig_pseudogene = ''
        pseudogene_length = sum(len(read) for read in self.reads)
        unitig_pseudogene_length = sum(len(unitig) for unitig in self.unitigs)
        self.position_to_read = np.zeros(pseudogene_length, dtype=np.int32)
        self.position_to_ug =np.zeros(unitig_pseudogene_length, dtype=np.int32)
        start = 0
        i = 0
        for idx,read in enumerate(self.reads):
        # for read in self.reads:
            # pseudogene += read
            end = start + len(read) - 1
            # self.sequence_info.append((start, end))
            self.seq_start.add(start)
            self.seq_end.add(end)
            # self.position_to_read[start:end+1] = idx
            start = end + 1
            i += 1
            # for _ in read:
            #     self.position_to_read.append(idx)
        self.pseudogene = ''.join(self.reads)
        # self.pseudogene = pseudogene
        i = 0
        start = 0
        for idx,unitig in enumerate(self.unitigs):
        # for unitig in self.unitigs:
            # unitig_pseudogene += unitig
            end = start + len(unitig) - 1
            # self.unitigs_info1.append((start, end))
            self.unitigs_info[idx] = start
            self.ug_start.add(start)
            self.ug_end.add(end)
            self.position_to_ug[start:end+1] = idx
            start = end + 1
            
            # self.unitigs_info[idx] = start
            # self.ug_start[start] = i
            # self.ug_end[end] = i
            # start = end + 1
            # for _ in unitig:
            #     self.position_to_ug.append(idx)
            i += 1
        unitig_pseudogene = ''.join(self.unitigs)
        seq_to_fasta([self.pseudogene],'seqfile.fa','seq')
        seq_to_fasta([unitig_pseudogene],'ugfile.fa','seq')


    def position_to_sequence(self, position):
        # 获取所有 unitig 结束位置的列表，用于二分查找
        ends = [end for _, end in self.sequence_info]
        # 二分查找定位到位置的区间
        index = bisect_right(ends, position)
        if index < len(self.sequence_info):
            start, _ = self.sequence_info[index]
            return index, start 
        return -1, -1
    
    def position_to_unitig(self, position):
        # 获取所有 unitig 结束位置的列表，用于二分查找
        # temp1 = time.time()
        # ends = [end+1 for _, end in self.unitigs_info1]
        # # 二分查找定位到位置的区间
        # index = bisect_right(ends, position)
        # if index < len(self.unitigs_info1):
        #     start, _ = self.unitigs_info1[index]
        #     temp2 = time.time()
        #     # self.t1 += temp2 - temp1
        #     return index, start
        # return -1, -1
        idx = self.position_to_ug[position]
        start = self.unitigs_info[idx]
        return idx,start
        
    
    def unitigs_suffix_prefix(self):
        # self.unitigs_prefix = {k: {} for k in range(self.min_len, self.kmer_len + 1)}
        # self.unitigs_suffix = {k: {} for k in range(self.min_len, self.kmer_len + 1)}
        # self.rc_unitigs_prefix = {k: {} for k in range(self.min_len, self.kmer_len + 1)}
        # self.rc_unitigs_suffix = {k: {} for k in range(self.min_len, self.kmer_len + 1)}

        for unitig in self.unitigs:
            rc_unitig = DNA_rev_complement(unitig)  # 计算反向互补序列

            # for k in range(self.min_len, self.kmer_len + 1):
            prefix = unitig[:self.kmer_len-1]
            suffix = unitig[-self.kmer_len+1:]
            # prefix = unitig[:self.kmer_len]
            # suffix = unitig[-self.kmer_len:]

            self.unitigs_prefix.add(prefix)
            self.unitigs_suffix.add(suffix)
            # self.unitigs_prefix.append(prefix)
            # self.unitigs_suffix.append(suffix)
            # self.rc_unitigs_prefix.append(rc_prefix)
            # self.rc_unitigs_suffix.append(rc_suffix)

    def unitigs_to_graph(self):
        for unitig in self.unitigs.keys():
            self.graph.add_node(unitig,coverage=self.unitigs[unitig])
        for unitig in self.unitigs.keys():
            prefix_kmer = unitig[:self.kmer_len-1]
            suffix_kmer = unitig[-(self.kmer_len-1):]
            if suffix_kmer in self.unitigs_prefix[self.kmer_len]:
                for unitig1 in self.unitigs_prefix[self.kmer_len][suffix_kmer]:
                    if unitig != unitig1:
                        self.graph.add_edge(unitig,unitig1,overlap=self.kmer_len-1)
            if prefix_kmer in self.unitigs_suffix[self.kmer_len]:
                for unitig1 in self.unitigs_suffix[self.kmer_len][prefix_kmer]:
                    if unitig != unitig1:
                        self.graph.add_edge(unitig1,unitig,overlap=self.kmer_len-1)

    def extend_unitigs1(self):
        original_nodes = list(self.graph.nodes())
        node_degrees = {node: (self.graph.in_degree(node), self.graph.out_degree(node)) for node in original_nodes}
        for node in original_nodes:
            in_degree, out_degree = node_degrees[node]
            if in_degree == 0 and out_degree == 0:
                self.extend_left(node)
                self.extend_right(node)
            if in_degree == 0:
                self.extend_left(node)
            if out_degree == 0:
                self.extend_right(node)
    
    def extend_unitigs(self):
        for unitig in self.unitigs:
            unitig_do = UnitigInfo()
            self.unitigs_ex_info[unitig] = unitig_do
        
        # self.parallel_unitig_extension(8)
        for unitig in self.unitigs:
            suffix = unitig[len(unitig)-self.kmer_len+1:]
            prefix = unitig[:self.kmer_len-1]
            rc_suffix = DNA_rev_complement(unitig[len(unitig)-self.kmer_len+1:])
            rc_prefix = DNA_rev_complement(unitig[:self.kmer_len-1]) 
            if suffix not in self.unitigs_prefix and rc_suffix not in self.unitigs_suffix:
                # if self.unitigs_dict[unitig] > 1:
                    self.change_order_right(unitig)
            if prefix not in self.unitigs_suffix and rc_prefix not in self.unitigs_prefix:
                # if self.unitigs_dict[unitig] > 1:
                    self.change_order_left(unitig)
        self.split_unitig()
        # self.add_extend_kmer()
    
    def add_extend_kmer(self):
        for unitig in self.unitigs_ex_info.keys():
            s_left = self.unitigs_ex_info[unitig].left_score
            s_right = self.unitigs_ex_info[unitig].right_score
            ex_left = self.unitigs_ex_info[unitig].left_extension
            ex_right = self.unitigs_ex_info[unitig].right_extension
            ex_left_list = self.unitigs_ex_info[unitig].left_extension_list
            ex_right_list = self.unitigs_ex_info[unitig].right_extension_list
            left_kmer = unitig[:self.kmer_len-1]
            right_kmer = unitig[len(unitig)-self.kmer_len+1:]
            for i in range(0,len(ex_left_list)):
                kmer = ex_left[len(ex_left_list)-i-1:len(ex_left_list)-i] + left_kmer
                left_kmer = kmer[:self.kmer_len-1]
                self.add_kmer(kmer,ex_left_list[i])
            for i in range(0,len(ex_right_list)):
                kmer = right_kmer + ex_right[i:i+1]
                right_kmer = kmer[:self.kmer_len-1]
                self.add_kmer(kmer,ex_right_list[i])

            

    def split_unitig(self):  
        for unitig in self.unitigs_ex_info.keys():
            s_left = self.unitigs_ex_info[unitig].left_score
            s_right = self.unitigs_ex_info[unitig].right_score
            ex_left = self.unitigs_ex_info[unitig].left_extension
            ex_right = self.unitigs_ex_info[unitig].right_extension
            
            prev = 0
            fragments = []
            fragments_cov = {}
            breakpoints = self.unitigs_ex_info[unitig].ug_break_point
            breakpoints.sort()
            flag = True
            for pos in breakpoints:
                if flag:
                    fragments_cov[ex_left+unitig[prev:pos+self.kmer_len-1]] = (s_left+len(unitig[prev:pos+self.kmer_len-1])*self.unitigs_dict[unitig])/(len(ex_left)+len(unitig[prev:pos+self.kmer_len-1]))
                else:
                    fragments_cov[unitig[prev:pos+self.kmer_len-1]] = self.unitigs_dict[unitig]
                fragments.append(unitig[prev:pos+self.kmer_len-1])
                prev = pos
                flag = False
            if flag:
                fragments_cov[ex_left+unitig+ex_right] = (s_left+s_right+len(unitig)*self.unitigs_dict[unitig])/(len(ex_left)+len(ex_right)+len(unitig))
            else:
                fragments_cov[unitig[prev:]+ex_right] = (s_right+len(unitig[prev:])*self.unitigs_dict[unitig])/(len(ex_right)+len(unitig[prev:]))
            self.new_unitigs.update(fragments_cov)
                
        
            
    def compute_fragments_cov(self,fragments):
        fragments_cov = {}
        for fragment in fragments:
            cov = 0
            for i in range(0,len(fragment)-self.kmer_len+1):

                res,rc_res = self.find_positions(fragment[i:i+self.kmer_len])
                cov += len(res)+len(rc_res)
            fragments_cov[fragment] = cov/(len(fragment)-self.kmer_len+1)
        return fragments_cov
    
    # def find_droplet_DNA(self, index, data_bytes_len):
    #     adrop = DNADroplet(bytes(data_bytes_len))
    #     adrop.head_index = index
    #     # adrop.get_tail_index()
    #     indexa = adrop.get_head_index_dna()
    #     # indexb = adrop.get_tail_index_dna()
    #     self.get_crc_path(indexa, data_bytes_len * 4 + adrop.crc_len, adrop.crc_len)

    # def get_crc_path(self, indexa, path_len, crc_len=8):
    #     self.crc_paths = []
    #     # if self.both_search:
    #     #     self.findSimPaths(indexa, indexb, path_len)
    #     # else:
    #     # temp1 = time.time()
    #     self.find_paths(indexa, path_len)
        
    #     for path in self.pathAB:
    #         pF_len = len(self.primerF)
    #         pathTrimPrimer = path[pF_len:]

    #         if self.crc_path_check(pathTrimPrimer, len(indexa), path_len, crc_len):
    #             # self.crc_paths.append(pathTrimPrimer)
    #             self.crc_paths.append(path+self.primerE)
    #         # else:
    #         #     print(len(self.pathAB))
    #         #     print("Wrong CRC Path detected")
    #     # temp3 = time.time()
    #     # self.t3 += temp3 - temp2

    # def crc_path_check(self, dnastr, indexa_dna_length, data_dna_length, crc_dna_length=8):
    #     # print('hello')
    #     #print(dnastr)
    #     # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
    #     # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
    #     # exit()
    #     data_bytes = DNAToBytes(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
    #     crc_bytes = DNAToBytes(
    #         dnastr[indexa_dna_length + data_dna_length - crc_dna_length:indexa_dna_length + data_dna_length])

    #     if crc16pure.crc16xmodem(data_bytes) == int.from_bytes(crc_bytes, byteorder='big', signed=False):
    #         return True
    #     else:
    #         return False
    
    def build_new_unitigs(self,unitig_file):
        # self.unitigs = fasta_to_dict(unitig_file)
        self.unitigs_dict = fasta_to_dict(unitig_file)
        self.unitigs = list(self.unitigs_dict.keys())
        self.update_prefix_suffix()
        self.unitigs_info = {}
        unitig_pseudogene_length = sum(len(unitig) for unitig in self.unitigs)
        self.position_to_ug =np.zeros(unitig_pseudogene_length, dtype=np.int32)
        # self.position_to_ug = []
        unitig_pseudogene = ''
        self.ug_start=set()
        self.ug_end=set()
        start = 0
        i = 0
        for idx,unitig in enumerate(self.unitigs):
            
            # unitig_pseudogene += unitig
            end = start + len(unitig) - 1
            # self.unitigs_info1.append((start, end))
            self.unitigs_info[idx] = start
            self.ug_start.add(start)
            self.ug_end.add(end)
            self.position_to_ug[start:end+1] = idx
            start = end + 1
            # self.unitigs_info[idx] = start
            # self.ug_start[start] = i
            # self.ug_end[end] = i
            # start = end + 1
            # for _ in unitig:
            #     self.position_to_ug.append(idx)
            i += 1
        unitig_pseudogene = ''.join(self.unitigs)
        seq_to_fasta([unitig_pseudogene],'newugfile.fa','seq')
        os.system('./kISS/build/kISS suffix_sort newugfile.fa -k 256 -t 4 --verbose')
        os.system('./kISS/build/kISS fmindex_build newugfile.fa -k 256 -t 4')
        self.ug_fmi = fm_index.FMIndex_Uint32_KISS1(self.ug_end)
        self.ug_fmi.load("/home/lzz/Projects/DBGPS_Python/newugfile.fa.fmi")
    
    def update_prefix_suffix(self):
        self.unitigs_prefix = {}
        self.unitigs_suffix = {}
        # for unitig in self.new_unitigs.keys():
        for unitig in self.unitigs:
            prefix = unitig[:self.kmer_len-1]
            suffix = unitig[-self.kmer_len+1:]
            if prefix not in self.unitigs_prefix:
                self.unitigs_prefix[prefix] = [unitig[self.kmer_len-1:]] 
            else:
                self.unitigs_prefix[prefix].append(unitig[self.kmer_len-1:])
            if suffix not in self.unitigs_suffix:
                self.unitigs_suffix[suffix] = [unitig[:len(unitig)-self.kmer_len+1]]
            else:
                self.unitigs_suffix[suffix].append(unitig[:len(unitig)-self.kmer_len+1])

    def merge_unitigs(self,unitig_file,output_file):
        # 存储已访问的 unitig 防止循环
        self.unitigs = fasta_to_dict(unitig_file)
        self.update_prefix_suffix()

        visited = set()
        
        # 最终合并的 contigs
        merged_contigs = {}
        new_unitigs_prefix = {}
        new_unitigs_suffix = {}
        merged_unitigs = []
        sorted_unitigs = sorted(self.unitigs, key=len, reverse=True)
        # 遍历所有 unitig
        
        for unitig in sorted_unitigs:
            if unitig in visited:
                continue  # 如果该 unitig 已经被合并，跳过

            visited.add(unitig)
            left_ex,left_count,left_coverage = self.left_merge(unitig,visited)
            right_ex,right_count,right_coverage = self.right_merge(unitig,visited)
                
            new_unitig = left_ex+unitig+right_ex
            merged_unitigs.append((new_unitig, len(new_unitig)))  # 保存 unitig 及其长度
        with open(output_file, 'w') as unitig_out:
            for i, (unitig, length) in enumerate(merged_unitigs, start=1):
                unitig_out.write(f">unitig{i} length={length}\n")
                unitig_out.write(f"{unitig}\n")
        #         prefix = new_unitig[:self.kmer_len-1]
        #         suffix = new_unitig[-self.kmer_len+1:]
        #         if prefix not in new_unitigs_prefix:
        #             new_unitigs_prefix[prefix] = [new_unitig[self.kmer_len-1:]]
        #         else:
        #             new_unitigs_prefix[prefix].append(new_unitig[self.kmer_len-1:])
        #         if suffix not in new_unitigs_suffix:
        #             new_unitigs_suffix[suffix] = [new_unitig[:len(new_unitig)-self.kmer_len+1]]
        #         else:
        #             new_unitigs_suffix[suffix].append(new_unitig[:len(new_unitig)-self.kmer_len+1])
        #         merged_contigs[new_unitig] = (left_coverage+self.new_unitigs[unitig]+right_coverage)/(left_count+right_count+1)
        # self.new_unitigs = merged_contigs
        # self.unitigs_prefix = new_unitigs_prefix
        # self.unitigs_suffix = new_unitigs_suffix
                
    # def merge_unitigs(self,unitig_file,output_file):
    #     # 存储已访问的 unitig 防止循环
    #     self.unitigs = fasta_to_dict(unitig_file)
    #     self.update_prefix_suffix()

    #     visited = set()
        
    #     # 最终合并的 contigs
    #     merged_contigs = {}
    #     new_unitigs_prefix = {}
    #     new_unitigs_suffix = {}
    #     # 遍历所有 unitig
        
    #     for unitig in self.new_unitigs.keys():
    #         if unitig in visited:
    #             continue  # 如果该 unitig 已经被合并，跳过

    #         visited.add(unitig)
    #         left_ex,left_count,left_coverage = self.left_merge(unitig,visited)
    #         right_ex,right_count,right_coverage = self.right_merge(unitig,visited)
                
    #         new_unitig = left_ex+unitig+right_ex
        
    #         prefix = new_unitig[:self.kmer_len-1]
    #         suffix = new_unitig[-self.kmer_len+1:]
    #         if prefix not in new_unitigs_prefix:
    #             new_unitigs_prefix[prefix] = [new_unitig[self.kmer_len-1:]]
    #         else:
    #             new_unitigs_prefix[prefix].append(new_unitig[self.kmer_len-1:])
    #         if suffix not in new_unitigs_suffix:
    #             new_unitigs_suffix[suffix] = [new_unitig[:len(new_unitig)-self.kmer_len+1]]
    #         else:
    #             new_unitigs_suffix[suffix].append(new_unitig[:len(new_unitig)-self.kmer_len+1])
    #         merged_contigs[new_unitig] = (left_coverage+self.new_unitigs[unitig]+right_coverage)/(left_count+right_count+1)
    #     self.new_unitigs = merged_contigs
    #     self.unitigs_prefix = new_unitigs_prefix
    #     self.unitigs_suffix = new_unitigs_suffix

    def left_merge(self, unitig, visited):
        merge_path = ''
        count = 0
        coverage = 0
        while True:

            prefix = unitig[:self.kmer_len-1]  # 当前 unitig 的后缀
            rc_prefix = DNA_rev_complement(prefix)
    #         if (len(self.unitigs_prefix.get(prefix, [])) + len(self.unitigs_suffix.get(rc_prefix, [])) == 1 and
    # len(self.unitigs_suffix.get(prefix, [])) + len(self.unitigs_prefix.get(rc_prefix, [])) == 1):
    #             if len(self.unitigs_suffix.get(prefix,[])) == 1:
            next_suffix_extends = self.unitigs_suffix.get(prefix, [])
            next_prefix_extends = self.unitigs_prefix.get(rc_prefix, [])
            max_contig = ''
            b_f = 0
            for next_suffix_extend in next_suffix_extends:
                next_unitig = next_suffix_extend+prefix
                if next_unitig not in visited and len(max_contig) < len(next_unitig):
                    b_f = 0
                    max_contig = next_suffix_extend
                    # merge_path = self.unitigs_suffix[prefix][0] + merge_path
                    # visited.add(self.unitigs_suffix[prefix][0] + prefix)
                    # unitig = self.unitigs_suffix[prefix][0] + prefix
                    # coverage += self.new_unitigs[unitig]
                    # count += 1
            for next_prefix_extend in next_prefix_extends:
                next_unitig = rc_prefix+next_prefix_extend
                if next_unitig not in visited and len(max_contig) < len(next_unitig):
                    b_f = 1
                    max_contig = next_prefix_extend
                # else:
                #     merge_path = DNA_rev_complement(self.unitigs_prefix[rc_prefix][0]) + merge_path
                #     visited.add(rc_prefix + self.unitigs_prefix[rc_prefix][0])
                #     coverage += self.new_unitigs[rc_prefix + self.unitigs_prefix[rc_prefix][0]] 
                #     unitig = DNA_rev_complement(rc_prefix + self.unitigs_prefix[rc_prefix][0])
                #     count += 1
            if max_contig == '':
                break
            if b_f==0:
                merge_path = max_contig + merge_path
                visited.add(max_contig+prefix)
                unitig = max_contig+prefix
                count += 1
            else:
                merge_path = DNA_rev_complement(max_contig)+merge_path
                visited.add(rc_prefix+max_contig)
                unitig = DNA_rev_complement(rc_prefix+max_contig)
                count += 1
            # else:
            #     break
        return merge_path,count,coverage
    
    # def left_merge(self, unitig, visited):
    #     merge_path = ''
    #     count = 0
    #     coverage = 0
    #     while True:

    #         prefix = unitig[:self.kmer_len-1]  # 当前 unitig 的后缀
    #         rc_prefix = DNA_rev_complement(prefix)
    #         if (len(self.unitigs_prefix.get(prefix, [])) + len(self.unitigs_suffix.get(rc_prefix, [])) == 1 and
    # len(self.unitigs_suffix.get(prefix, [])) + len(self.unitigs_prefix.get(rc_prefix, [])) == 1):
    #             if len(self.unitigs_suffix.get(prefix,[])) == 1:
    #                 merge_path = self.unitigs_suffix[prefix][0] + merge_path
    #                 visited.add(self.unitigs_suffix[prefix][0] + prefix)
    #                 unitig = self.unitigs_suffix[prefix][0] + prefix
    #                 coverage += self.new_unitigs[unitig]
    #                 count += 1
    #             else:
    #                 merge_path = DNA_rev_complement(self.unitigs_prefix[rc_prefix][0]) + merge_path
    #                 visited.add(rc_prefix + self.unitigs_prefix[rc_prefix][0])
    #                 coverage += self.new_unitigs[rc_prefix + self.unitigs_prefix[rc_prefix][0]] 
    #                 unitig = DNA_rev_complement(rc_prefix + self.unitigs_prefix[rc_prefix][0])
    #                 count += 1
    #         else:
    #             break
    #     return merge_path,count,coverage
    
    def right_merge(self, unitig, visited):
        merge_path = ''
        count = 0
        coverage = 0
        while True:

            suffix = unitig[-self.kmer_len+1:]  # 当前 unitig 的后缀
            rc_suffix = DNA_rev_complement(suffix)
    #         if (len(self.unitigs_suffix.get(suffix, [])) + len(self.unitigs_prefix.get(rc_suffix, [])) == 1 and
    # len(self.unitigs_prefix.get(suffix, [])) + len(self.unitigs_suffix.get(rc_suffix, [])) == 1):
    #             if len(self.unitigs_prefix.get(suffix,[])) == 1:
            next_prefix_extends = self.unitigs_prefix.get(suffix, [])
            next_suffix_extends = self.unitigs_suffix.get(rc_suffix, [])
            max_contig = ''
            b_f = 0
            for next_prefix_extend in next_prefix_extends:
                next_unitig = suffix+next_prefix_extend
                if next_unitig not in visited and len(max_contig) < len(next_unitig):
                    b_f = 0
                    max_contig = next_prefix_extend
                    # merge_path = merge_path + self.unitigs_prefix[suffix][0]
                    # visited.add(suffix + self.unitigs_prefix[suffix][0])
                    # unitig = suffix + self.unitigs_prefix[suffix][0]
                    # coverage += self.new_unitigs[unitig]
                    # count += 1
            for next_suffix_extend in next_suffix_extends:
                next_unitig = next_suffix_extend+rc_suffix
                if next_unitig not in visited and len(max_contig) < len(next_unitig):
                    b_f = 1
                    max_contig = next_suffix_extend
                # else:
                #     merge_path = merge_path + DNA_rev_complement(self.unitigs_suffix[rc_suffix][0])
                #     visited.add(self.unitigs_suffix[rc_suffix][0] + rc_suffix)
                #     coverage += self.new_unitigs[self.unitigs_suffix[rc_suffix][0] + rc_suffix]
                #     unitig = DNA_rev_complement(self.unitigs_suffix[rc_suffix][0] + rc_suffix)
                #     count += 1
            if max_contig == '':
                break
            if b_f==0:
                merge_path = merge_path+max_contig
                visited.add(suffix+max_contig)
                unitig = suffix+max_contig
                count += 1
            else:
                merge_path = merge_path+DNA_rev_complement(max_contig)
                visited.add(max_contig+rc_suffix)
                unitig = DNA_rev_complement(max_contig+rc_suffix)
                count += 1
            # else:
            #     break
        return merge_path,count,coverage


    def find_paths(self,indexa, path_len):
        self.pathsA = []
        self.pathAB = []
        self.pathABLen = len(indexa) + path_len
        res,rc_res = self.find_ug_positions((self.primerF + indexa)[-self.kmer_len+1:])
        # res = self.del_unitig_linkpoint(self.kmer_len-1,res)
        # rc_res = self.del_unitig_linkpoint(self.kmer_len-1,rc_res)
        for pos in res:
            ug_idx,start = self.position_to_unitig(pos)
            if pos - start == len(self.unitigs[ug_idx])-self.kmer_len+1:
                continue
            if len(self.unitigs[ug_idx][pos-start+self.kmer_len-1:]) >= path_len:
                self.pathAB.append((self.primerF + indexa + self.unitigs[ug_idx][pos-start+self.kmer_len-1:])[:self.pathABLen + len(self.primerF)])
            else:
                self.pathsA.append(self.primerF + indexa + self.unitigs[ug_idx][pos-start+self.kmer_len-1:])
        for pos in rc_res:
            ug_idx,start = self.position_to_unitig(pos)
            if pos == start:
                continue
            if len(self.unitigs[ug_idx][:pos-start]) >= path_len:
                self.pathAB.append((self.primerF + indexa + DNA_rev_complement(self.unitigs[ug_idx][:pos-start]))[:self.pathABLen + len(self.primerF)])
            else:
                self.pathsA.append(self.primerF + indexa + DNA_rev_complement(self.unitigs[ug_idx][:pos-start]))
        
        while len(self.pathAB) < 100:
            
            tmpPaths = []
            for path in self.pathsA:
                kmer = path[-self.kmer_len+1:]
                ex = self.find_extend(kmer)
                for temp in ex:
                    if len(temp) == 0:
                        continue
                    if len(temp) + len(path)>=self.pathABLen + len(self.primerF):
                        self.pathAB.append((path+temp)[:self.pathABLen + len(self.primerF)])
                    else:
                        tmpPaths.append(path+temp)
                # res,rc_res = self.find_ug_positions(kmer)
                # res = self.del_linkpoint_pos1(self.kmer_len-1,res)
                # for pos in res:
                #     ex = ''
                #     if len(ex) == 0:
                #         continue
                #     if len(ex) + len(path)>=self.pathABLen + len(self.primerF):
                #         self.pathAB.append((path+ex)[:self.pathABLen + len(self.primerF)])
                #     else:
                #         tmpPaths.append(path+ex)
                # rc_res = self.del_linkpoint_pos1(self.kmer_len-1,rc_res)
                # for pos in rc_res:
                #     ex = self.find_backward(pos)
                #     if len(ex) == 0:
                #         continue
                #     if len(ex) + len(path)>=self.pathABLen + len(self.primerF):
                #         self.pathAB.append((path+DNA_rev_complement(ex))[:self.pathABLen + len(self.primerF)])
                #     else:
                #         tmpPaths.append(path+DNA_rev_complement(ex))
            
            if len(tmpPaths) > 0:
                self.pathsA = tmpPaths
            else:
                return
    
    def find_extend(self,kmer):
        ex = []
        rc_kmer = DNA_rev_complement(kmer)
        if kmer in self.unitigs_prefix:
            ex = self.unitigs_prefix[kmer]
        if rc_kmer in self.unitigs_suffix:
            ex = []
            temps = []
            temps = self.unitigs_suffix[rc_kmer]
            for temp in temps:
                ex.append(DNA_rev_complement(temp))
        return ex

    def find_forward(self,pos):
        ug_idx,start = self.position_to_unitig(pos)
        if pos == start:
            print(f"unitig id {ug_idx},extend {self.unitigs[ug_idx][self.kmer_len-1:]}")
            return self.unitigs[ug_idx][self.kmer_len-1:]
        return ''
    
    def find_backward(self,pos):
        ug_idx,start = self.position_to_unitig(pos)
        if pos - start == len(self.unitigs[ug_idx])-self.kmer_len+1:
            print(f"unitig id {ug_idx},extend {self.unitigs[ug_idx][:len(self.unitigs[ug_idx])-self.kmer_len+1]}")
            return self.unitigs[ug_idx][:len(self.unitigs[ug_idx])-self.kmer_len+1]
        return ''

    def extend_left(self,node):
        # temp1 = time.time()
        currK = self.kmer_len
        self.currNodes = [node]
        i = 0
        while i<=self.dataLen and len(self.currNodes) < self.max_path_num:
            tmpNodes = self.change_order_left(currK)
            if len(tmpNodes)==0:
                break
            self.currNodes = tmpNodes
            i+=1
        # temp2 = time.time()
        # self.findright1 += temp2-temp1 
            
    def extend_right(self,node):
        # temp1 = time.time()
        currK = self.kmer_len
        self.currNodes = [node]
        i = 0
        while i<=self.dataLen and len(self.currNodes) < self.max_path_num:
            tmpNodes = self.change_order_right(currK)
            if len(tmpNodes)==0:
                break
            self.currNodes = tmpNodes
            i+=1
        # temp2 = time.time()
        # self.findright1 += temp2-temp1 

    
    def extend_node(self, indices, positions, node, threshold, isleft):
        result = {}
        paths = {}
        paths[node]=len(indices)
        path_indices = {node: list(range(len(indices)))}
        for _ in range(len(node), self.kmer_len):
            tmpPaths = {}
            for path in paths:
                curr_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
                bases_above_threshold = self.score_bases(curr_bases,indices,positions,threshold,path_indices[path],isleft)
                if len(bases_above_threshold) == 0:
                    result[path] = len(path_indices[path])
                    continue
                del path_indices[path]
                for base in bases_above_threshold.keys():
                    if not isleft:
                        tmpPaths[path+base]=len(bases_above_threshold[base])
                        path_indices[path+base] = bases_above_threshold[base]
                    else:
                        tmpPaths[base+path]=len(bases_above_threshold[base])
                        path_indices[base+path] = bases_above_threshold[base]
            if len(tmpPaths) > 0:
                paths = tmpPaths
            else:
                return result
        result.update(paths)
        return result


    def r_record_break_point(self,fb,kmer_position):
        ug,start = self.position_to_unitig(kmer_position)
        if fb == 0:
            # if kmer_position != start+len(self.unitigs[ug])-self.kmer_len:
            if kmer_position != start:
                if kmer_position-start not in self.unitigs_ex_info[self.unitigs[ug]].ug_break_point:
                    self.unitigs_ex_info[self.unitigs[ug]].ug_break_point.append(kmer_position-start)
        else:
            # if kmer_position != start:
            if kmer_position != start+len(self.unitigs[ug])-self.kmer_len:
                if kmer_position-start+1 not in self.unitigs_ex_info[self.unitigs[ug]].ug_break_point:
                    self.unitigs_ex_info[self.unitigs[ug]].ug_break_point.append(kmer_position-start+1)

    def l_record_break_point(self,fb,kmer_position):
        ug,start = self.position_to_unitig(kmer_position)
        if fb == 0:
            # if kmer_position != start:
            if kmer_position != start+len(self.unitigs[ug])-self.kmer_len:
                if kmer_position-start+1 not in self.unitigs_ex_info[self.unitigs[ug]].ug_break_point:
                    self.unitigs_ex_info[self.unitigs[ug]].ug_break_point.append(kmer_position-start+1)
        else:
            # if kmer_position != start+len(self.unitigs[ug])-self.kmer_len:
            if kmer_position != start:
                if kmer_position-start not in self.unitigs_ex_info[self.unitigs[ug]].ug_break_point:
                    self.unitigs_ex_info[self.unitigs[ug]].ug_break_point.append(kmer_position-start)

        # res,rc_res = self.find_ug_positions(kmer)
        # res,rc_res = self.del_linkpoint_pos(len(kmer),res,rc_res)
        # res = res + rc_res
        # result = 0
        # ug_idx = ''
        # for pos in res:
        #     ug,start = self.position_to_unitig(pos)
        #     if pos + len(kmer) < start+len(self.unitigs[ug])-1 and pos > start:
        #         result = pos
        #         ug_idx = ug
        # return result,ug_idx
                
        # if len(res)>1:
        #     return
        # if res[0] not in self.unitigs_ex_info[self.unitigs[ug_idx]].ug_break_point:
        #     self.unitigs_ex_info[self.unitigs[ug_idx]].ug_break_point.append(res[0])
        # return res,ug_idx
            
    def del_linkpoint_pos(self,k,res,rc_res):
        new_res = []
        new_rc_res = []
        for pos in res:
            ug,start = self.position_to_ug(pos)
            if pos + k <= start + len(self.unitigs[ug]):
                new_res.append(pos)
        for pos in rc_res:
            ug,start = self.position_to_ug(pos)
            if pos + k <= start + len(self.unitigs[ug]):
                new_rc_res.append(pos)
        return new_res,new_rc_res
    
    def del_seq_linkpoint(self,k,res):
        new_res = []
        for pos in res:
            end = pos + k - 1
            idx1 = self.position_to_read[pos]
            idx2 = self.position_to_read[end]
            if idx1 != idx2:
                continue
            new_res.append(pos)
            # flag = True
            # for i in range(0,k-1):
            #     if pos+i in self.seq_end:
            #         flag = False
            # if flag:
            #     new_res.append(pos)

        return new_res
    
    def del_unitig_linkpoint(self,k,res):
        new_res = []
        for pos in res:
            end = pos + k - 1
            if end >= len(self.position_to_ug):
                continue
            idx1 = self.position_to_ug[pos]
            idx2 = self.position_to_ug[end]
            if idx1 != idx2:
                continue
            new_res.append(pos)
            # flag = True
            # for i in range(0,k-1):
            #     if pos+i in self.ug_end:
            #         flag = False
            # if flag:
            #     new_res.append(pos)
        return new_res


    def find_positions(self,node):
        
        res = []
        rc_res = []
        result = self.fmi.get_range(node,0)
        curr_offset = abs(result[0]-result[1])
        rc_result = self.fmi.get_range(DNA_rev_complement(node),0)
        rc_curr_offset = abs(rc_result[0]-rc_result[1])
        temp1 = time.time()
        if rc_curr_offset+curr_offset>=self.mid_threshold:
            res = self.fmi.get_offsets(result[0],result[1])
            rc_res = self.fmi.get_offsets(rc_result[0],rc_result[1])
        # temp2 = time.time()
        # self.t3 += temp2 - temp1
        return res,rc_res
    
    def find_ug_positions(self,node):
        temp1 = time.time()
        # self.t3 += temp2 - temp1
        result = self.ug_fmi.get_range(node,0)
        res = self.ug_fmi.get_offsets(result[0],result[1])
        rc_result = self.ug_fmi.get_range(DNA_rev_complement(node),0)
        rc_res = self.ug_fmi.get_offsets(rc_result[0],rc_result[1])
        temp2 = time.time()
        self.t1 += temp2 - temp1
        return res,rc_res
    
    def check_tail_head(self,kmer,rc_kmer):
        if kmer not in self.unitigs_prefix and kmer not in self.rc_unitigs_prefix and rc_kmer not in self.unitigs_suffix and rc_kmer not in self.rc_unitigs_suffix:
            return False
        else:
            return True
        
    def get_kmers(self, sequence, k):
        temp1 = time.time()
        kmers = Counter()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmers[kmer] += 1
        temp2 = time.time()
        # self.t3 += temp2-temp1
        return kmers

    # def get_kmers(self,sequence, k):
    #     # 使用位操作将序列转换为整数表示
    #     kmers = Counter()
    #     base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    #     kmer_value = 0
    #     for i in range(k):
    #         kmer_value = (kmer_value << 2) | base_map[sequence[i]]
    #     kmers[kmer_value] += 1
    #     for i in range(k, len(sequence)):
    #         kmer_value = ((kmer_value << 2) & ((1 << (2 * k)) - 1)) | base_map[sequence[i]]
    #         kmers[kmer_value] += 1
    #     return kmers

    def align_and_filter(self, pattern, texts, k, threshold):
        # self.t1 += 1
        temp1 = time.time()
        command = ['./WFA2-lib/examples/bin/wfa_bindings',str(k),str(threshold),pattern] + texts
        # command_str = ' '.join(command)
        result = subprocess.run(command, capture_output=True, text=True)
        temp2 = time.time()
        # self.t3 += temp2-temp1
        # 检查是否运行成功
        if result.returncode == 0:
            # 解析输出的比对分数，输出会是以空格分隔的分数
            valid_indices = list(map(int, result.stdout.strip().split()))
            return valid_indices
        else:
            print(f"Error occurred: {result.stderr}")
            return None

    def dna_alignment_score1(self, seq1, seq2):
        temp1 = time.time()
        m, n = len(seq1), len(seq2)
        
        # Initialize the DP table with dimensions (m+1) x (n+1)
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        
        # Fill the DP table
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if seq1[i - 1] == seq2[j - 1]:
                    # If characters match, take the diagonal value (no penalty for match)
                    dp[i][j] = dp[i - 1][j - 1]  # Keep the score same as previous
                else:
                    # If they do not match, consider insertion, deletion, or mismatch (all penalized by 1)
                    dp[i][j] = max(dp[i - 1][j] - 1,  # Deletion
                                dp[i][j - 1] - 1,  # Insertion
                                dp[i - 1][j - 1] - 1)  # Mismatch
        
        # The alignment score is the value in the bottom-right corner of the table
        temp2 = time.time()
        # self.t3 += temp2 - temp1
        return dp[m][n]

    def dna_alignment_score2(self, seq1, seq2):
        temp1 = time.time()
        self.wfa_aligner.alignEnd2End(seq1,seq2)
        result = self.wfa_aligner.getAlignmentScore()
        # result = self.wfa_aligner.jaccard_wfa(seq1,seq2,self.check_k,self.max_mis_num)
        temp2 = time.time()
        # self.t3 += temp2 - temp1
        return result

    def prefix_to_base(self,node,positions,currK):
        sum_num = 0
        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        temp1 = time.time()
        res,rc_res = self.find_positions(node[-currK+1:])
        temp2 = time.time()
        self.t3 += temp2 - temp1
        # res = self.del_seq_linkpoint(self.min_k-1,res)
        # rc_res = self.del_seq_linkpoint(self.min_k-1,rc_res)
        
        # res,rc_res = self.check_extend(res,rc_res)
        template_kmers = self.get_kmers(node[len(node)-self.kmer_len:len(node)-currK+1], self.check_k)
        for pos in res:
            # if pos + self.min_k-1 >= len(self.position_to_read) or pos-self.kmer_len+self.min_k-1 < 0 or self.position_to_read[pos-self.kmer_len+self.min_k-2] != self.position_to_read[pos + self.min_k-1]:
            if pos + currK-1 >= len(self.pseudogene) or pos-self.kmer_len+currK-1 < 0 or pos +currK-2 in self.seq_end:
                continue
            temp1 = time.time()
            seq_kmers = self.get_kmers(self.pseudogene[pos-self.kmer_len+currK-1:pos], self.check_k)
            # texts.append(self.pseudogene[pos-self.kmer_len+self.min_k-1:pos])
            # new_res.append(pos)
            intersection = sum((template_kmers & seq_kmers).values())
            jaccard_idx = intersection / (self.kmer_len-currK+1-self.check_k+1)
            temp2 = time.time()
            # self.t3 += temp2 - temp1
            score = self.max_mis_num +1
            if jaccard_idx >= 1-self.check_k*self.max_mis_num/(self.kmer_len-currK+1-self.check_k+1):
                score = -self.dna_alignment_score2(node[len(node)-self.kmer_len:len(node)-currK+1],self.pseudogene[pos-self.kmer_len+currK-1:pos])
            
            if score <= self.max_mis_num:
        # valid_indices = self.align_and_filter(pattern,texts,self.check_k,self.max_mis_num)
        # for idx in valid_indices:
        #         pos = new_res[idx]
                if self.pseudogene[pos+currK-1] != 'N':
                    key = self.pseudogene[pos + currK-1]
                    edgeScore[key] += 1
                    sum_num += 1
                    if sum_num > 10:
                        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        return edgeScore
                    if key not in positions:
                        positions[key] = [pos+currK-1]
                    else:
                        positions[key].append(pos+currK-1)
        texts = []
        for pos in rc_res:
            # if pos + self.kmer_len >= len(self.position_to_read) or pos-1 < 0 or self.position_to_read[pos+self.kmer_len-1] != self.position_to_read[pos-1]:
            if pos + self.kmer_len >= len(self.pseudogene) or pos-1 < 0 or pos in self.seq_start:
                continue
            temp1 = time.time()
            seq_kmers = self.get_kmers(DNA_rev_complement(self.pseudogene[pos+currK-1:pos+self.kmer_len]), self.check_k)
            # texts.append(DNA_rev_complement(self.pseudogene[pos+self.min_k-1:pos+self.kmer_len]))
            # new_rc_res.append(pos)
            intersection = sum((template_kmers & seq_kmers).values())
            jaccard_idx = intersection / (self.kmer_len-currK+1-self.check_k+1)
            temp2 = time.time()
            # self.t3 += temp2 - temp1
            score = self.max_mis_num +1
            if jaccard_idx >= 1-self.check_k*self.max_mis_num/(self.kmer_len-currK+1-self.check_k+1):
                score = -self.dna_alignment_score2(node[len(node)-self.kmer_len:len(node)-currK+1],DNA_rev_complement(self.pseudogene[pos+currK-1:pos+self.kmer_len]))

            if score <= self.max_mis_num:
        # valid_indices = self.align_and_filter(pattern,texts,self.check_k,self.max_mis_num)
        # for idx in valid_indices:
        #         pos = new_rc_res[idx]
                if self.pseudogene[pos-1] != 'N':
                    key = complement_base(self.pseudogene[pos-1])
                    edgeScore[key] += 1
                    sum_num += 1
                    if sum_num > 10:
                        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        return edgeScore
                    if key not in positions:
                        positions[key] = [pos-1]
                    else:
                        positions[key].append(pos-1)
        temp2 = time.time()
        # self.t1 += temp2 - temp1
        return edgeScore

    def suffix_to_base(self,node,positions,currK):
        
        sum_num = 0
        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        temp1 = time.time()
        res,rc_res = self.find_positions(node[:currK-1])
        temp2 = time.time()
        self.t3 += temp2 - temp1
        # res = self.del_seq_linkpoint(self.min_k-1,res)
        # rc_res = self.del_seq_linkpoint(self.min_k-1,rc_res)
        # res,rc_res = self.check_extend(res,rc_res)
        template_kmers = self.get_kmers(node[currK-1:self.kmer_len], self.check_k)

        for pos in res:

            # if pos + self.kmer_len >= len(self.position_to_read) or pos-1 < 0 or self.position_to_read[pos+self.kmer_len-1] != self.position_to_read[pos-1]:
            if pos + self.kmer_len >= len(self.pseudogene) or pos-1 < 0 or pos in self.seq_start:
                continue
            temp1 = time.time()
            seq_kmers = self.get_kmers(self.pseudogene[pos+currK-1:pos+self.kmer_len], self.check_k)
            # texts.append(self.pseudogene[pos+self.min_k-1:pos+self.kmer_len])
            # new_res.append(pos)
            intersection = sum((template_kmers & seq_kmers).values())
            jaccard_idx = intersection / (self.kmer_len-currK+1-self.check_k+1)
            temp2 = time.time()
            # self.t3 += temp2 - temp1
            score = self.max_mis_num +1
            if jaccard_idx >= 1-self.check_k*self.max_mis_num/(self.kmer_len-currK+1-self.check_k+1):
                score = -self.dna_alignment_score2(node[currK-1:self.kmer_len],self.pseudogene[pos+currK-1:pos+self.kmer_len])
            
            if score <= self.max_mis_num:
        # valid_indices = self.align_and_filter(pattern,texts,self.check_k,self.max_mis_num)
        # for idx in valid_indices:
        #         pos = new_res[idx]
                if self.pseudogene[pos-1] != 'N':
                    key = self.pseudogene[pos-1]
                    edgeScore[key] += 1
                    sum_num += 1
                    if sum_num > 10:
                        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        return edgeScore
                    if key not in positions:
                        positions[key] = [pos-1]
                    else:
                        positions[key].append(pos-1)
        texts = []
        for pos in rc_res:

            # if pos + self.min_k-1 >= len(self.position_to_read) or pos-self.kmer_len+self.min_k-1 < 0 or self.position_to_read[pos-self.kmer_len+self.min_k-2] != self.position_to_read[pos + self.min_k-1]:
            if pos + currK-1 >= len(self.pseudogene) or pos-self.kmer_len+currK-1 < 0 or pos +currK-2 in self.seq_end:
                continue
            temp1 = time.time()
            seq_kmers = self.get_kmers(DNA_rev_complement(self.pseudogene[pos-self.kmer_len+currK-1:pos]), self.check_k)
            # texts.append(DNA_rev_complement(self.pseudogene[pos-self.kmer_len+self.min_k-1:pos]))
            # new_rc_res.append(pos)
            intersection = sum((template_kmers & seq_kmers).values())
            jaccard_idx = intersection / (self.kmer_len-currK+1-self.check_k+1)
            
            temp2 = time.time()
            # self.t3 += temp2 - temp1
            score = self.max_mis_num +1
            if jaccard_idx >= 1-self.check_k*self.max_mis_num/(self.kmer_len-currK+1-self.check_k+1):
                score = -self.dna_alignment_score2(node[currK-1:self.kmer_len],DNA_rev_complement(self.pseudogene[pos-self.kmer_len+currK-1:pos]))

            if score <= self.max_mis_num:
        # valid_indices = self.align_and_filter(pattern,texts,self.check_k,self.max_mis_num)
        # for idx in valid_indices:
        #         pos = new_rc_res[idx]
                if self.pseudogene[pos+currK-1] != 'N':
                    key = complement_base(self.pseudogene[pos + currK-1])
                    edgeScore[key] += 1
                    sum_num += 1
                    if sum_num > 10:
                        edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        return edgeScore
                    if key not in positions:
                        positions[key] = [pos+currK-1]
                    else:
                        positions[key].append(pos+currK-1)
        temp2 = time.time()
        # self.t1 += temp2 - temp1
        return edgeScore

    def check_extend(self, res, rc_res):
        new_res = []
        new_rc_res = []
        for pos in res:
            # if pos +self.min_k-2 in self.seq_end:
            if pos +self.min_k-1 >= len(self.position_to_read) or self.position_to_read[pos +self.min_k-1] != self.position_to_read[pos +self.min_k-2]:
                continue
            new_res.append(pos)
        for pos in rc_res:
            # if pos in self.seq_start:
            if pos <= 0 or self.position_to_read[pos] != self.position_to_read[pos-1]:
                continue
            new_rc_res.append(pos)
        return new_res,new_rc_res

    def change_order_right(self,unitig):
        pathsA = {}
        paths_score_list = {}
        paths_score_list[unitig] = []
        pathsA[unitig] = [0,self.kmer_len]
        max_confidence_path = ''
        max_confidence = 0
        max_count = 0
        break_paths = {}
        max_break_info = {}
        i = 0
        while i < 100 and len(pathsA) < 50:
            tmpPaths = {}
            for node in pathsA:
                flag = False
                score_list = paths_score_list[node]
                currK = pathsA[node][1]
                prev_pos = []
                if len(pathsA[node]) > 2:
                    prev_pos = pathsA[node][2:]
                break_info = {}
                prev_offest = 0
                while currK >= self.min_k:
                    positions = {}
                    edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                    # currK = self.kmer_len if len(node) > self.kmer_len else len(node)
                    pNode = node[len(node)-currK:]
                    temp1 = time.time()
        
                    res,rc_res = self.find_positions(pNode[1:])
                    temp2 = time.time()
                    # self.t3 += temp2 - temp1
                    # res = self.del_seq_linkpoint(len(pNode)-1,res)
                    # rc_res = self.del_seq_linkpoint(len(pNode)-1,rc_res)
                    # print(res)
                    # print(rc_res)
                    for pos in res:
                        # if pos + currK-1 == len(self.position_to_read) or self.position_to_read[pos + currK-2] != self.position_to_read[pos + currK-1] or (currK < self.mid_k and pos+ currK-2 not in prev_pos):
                        if pos + currK-1 >= len(self.position_to_read) or pos + currK-2 in self.seq_end or (currK < self.mid_k and pos+ currK-2 not in prev_pos):
                            continue
                        if self.pseudogene[pos+currK-1] != 'N':
                            key = self.pseudogene[pos + currK-1]
                            edgeScore[key] += 1
                            if key not in positions:
                                positions[key] = [pos+currK-1]
                            else:
                                positions[key].append(pos+currK-1)
                    for pos in rc_res:
                        # if pos-1 == 0 or self.position_to_read[pos] != self.position_to_read[pos-1] or (currK < self.mid_k and pos not in prev_pos):
                        if pos in self.seq_start or (currK < self.mid_k and pos not in prev_pos):
                            continue
                        if self.pseudogene[pos-1] != 'N':
                            key = complement_base(self.pseudogene[pos-1])
                            edgeScore[key] += 1
                            if key not in positions:
                                positions[key] = [pos-1]
                            else:
                                positions[key].append(pos-1)
                    has_value_exceeding = any(value >= self.mid_threshold for value in edgeScore.values())
                    if currK <= self.mid_k and not has_value_exceeding:
                        # positions = {}
                        if currK > self.min_k:
                            if currK-self.gap_len >= self.min_k:
                                edgeScore = self.prefix_to_base(node,positions,currK-self.gap_len)
                            else:
                                edgeScore = self.prefix_to_base(node,positions,self.min_k)
                        currK = self.min_k
                    
                    for base, count in edgeScore.items():
                        if (count >= self.mid_threshold and currK >= self.mid_k) or count >= self.min_threshold:
                            if node in paths_score_list and max_confidence_path != node:
                                    del paths_score_list[node]
                            paths_score_list[node+base] = score_list+[(currK/self.kmer_len)*count]
                            if currK < self.kmer_len:
                                tmpPaths[node+base] = [pathsA[node][0]+(currK/self.kmer_len)*count,currK+1]
                                tmpPaths[node+base] = tmpPaths[node+base]+positions[base]
                            else:
                                kmer = node[len(node)-self.kmer_len+1:] + base
                                res,rc_res = self.find_ug_positions(kmer)
                                    # res = self.del_unitig_linkpoint(len(kmer),res)
                                    # rc_res = self.del_unitig_linkpoint(len(kmer),rc_res)
                                    # if len(res) + len(rc_res) == 0 or count < self.min_threshold:
                                if len(res) + len(rc_res) == 0:
                                    tmpPaths[node+base] = [pathsA[node][0]+(currK/self.kmer_len)*count,currK]
                                    tmpPaths[node+base] = tmpPaths[node+base]+positions[base]
                                        # score_list = paths_score_list[node]
                                else:
                                    flag = True
                                    # if kmer[1:] not in self.unitigs_suffix and DNA_rev_complement(kmer[1:]) not in self.unitigs_prefix:
                                    if len(res) != 0:
                                        break_info[base] = (0,res[0])
                                    else:
                                        break_info[base] = (1,rc_res[0])
                                    confidence = self.compute_path_confidence(node,pathsA[node][0],unitig)
                                    if max_confidence == 0:
                                        max_confidence_path = node
                                        max_count = pathsA[node][0]
                                        max_confidence = confidence
                                        max_break_info = break_info
                                    else:
                                        if confidence > max_confidence:
                                            max_confidence_path = node
                                            max_count = pathsA[node][0]
                                            max_confidence = confidence
                                            max_break_info = break_info
                    if len(tmpPaths) > 0 or len(break_info) > 0:
                        break
                    else:
                        currK -= 1
                # if len(break_info) > 0:
                #     confidence = self.compute_path_confidence(node,pathsA[node][0],unitig)
                #     if max_confidence == 0:
                #         max_confidence_path = node
                #         max_count = pathsA[node][0]
                #         max_confidence = confidence
                #         max_break_info = break_info
                #     else:
                #         if confidence > max_confidence:
                #             max_confidence_path = node
                #             max_count = pathsA[node][0]
                #             max_confidence = confidence
                #             max_break_info = break_info
                if currK < self.min_k:
                    break_paths[node] = (pathsA[node][0],self.kmer_len)
            if len(tmpPaths) > 0:
                pathsA = tmpPaths
            else:
                best_path = ''
                count = 0
                if len(max_confidence_path) > 0:
                    best_path = max_confidence_path
                    count = max_count
                    for base in max_break_info.keys():
                        kmer = best_path[len(best_path)-self.kmer_len+1:] + base
                        (fb,kmer_position) = max_break_info[base]
                        self.r_record_break_point(fb,kmer_position)
                    # self.unitigs_suffix.add(best_path[-self.kmer_len+1:])
                    # self.rc_unitigs_prefix.append(DNA_rev_complement(best_path[-self.kmer_len+1:]))
                else:
                    best_path,count = self.choose_best_path(break_paths,unitig)
                self.unitigs_suffix.discard(unitig[-self.kmer_len+1:])
                self.unitigs_suffix.add(best_path[-self.kmer_len+1:])
                self.unitigs_ex_info[unitig].right_extension = best_path[len(unitig):]
                self.unitigs_ex_info[unitig].right_score = count
                # temp = paths_score_list[best_path]
                # self.unitigs_ex_info[unitig].right_extension_list = paths_score_list[best_path]
                # for base in max_last_base:
                return best_path[-self.kmer_len+1:]
            i += 1
        break_paths.update(pathsA)
        best_path,count = self.choose_best_path(break_paths,unitig)
        self.unitigs_ex_info[unitig].right_extension = best_path[len(unitig):]
        self.unitigs_ex_info[unitig].right_score = count
        # self.unitigs_ex_info[unitig].right_extension_list = paths_score_list[best_path]
        return best_path[-self.kmer_len+1:]
        
        
            
    def change_order_left(self,unitig):
        pathsA = {}
        paths_score_list = {}
        paths_score_list[unitig] = []
        pathsA[unitig] = [0,self.kmer_len]
        max_confidence_path = ''
        max_confidence = 0
        max_count = 0
        break_paths = {}
        max_break_info = {}
        i = 0
        while i <100 and len(pathsA) < 50:
            tmpPaths = {}
            for node in pathsA:
                flag = False
                score_list = paths_score_list[node]
                currK = pathsA[node][1]
                break_info = {}
                prev_pos = []
                if len(pathsA[node]) > 2:
                    prev_pos = pathsA[node][2:]
                prev_offest = 0
                while currK >= self.min_k:
                    temp1 = time.time()
                    positions = {}
                    edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                    # currK = self.kmer_len if len(node) > self.kmer_len else len(node)
                    pNode = node[0:currK]
                    temp1 = time.time()
                    res,rc_res = self.find_positions(pNode[:len(pNode)-1])
                    temp2 = time.time()
                    # self.t3 += temp2 - temp1
                    # res = self.del_seq_linkpoint(len(pNode)-1,res)
                    # rc_res = self.del_seq_linkpoint(len(pNode)-1,rc_res)
                    # if currK == self.mid_k:
                    #     rc_res,res = self.check_extend(rc_res,res)
                    #     self.suffix_to_base(node,res,rc_res,positions,rc_positions,edgeScore)
                    # else:
                    for pos in res:
                            # check_pos = True
                            # for ref_pos in prev_pos:
                            #     if ref_pos-self.max_gap_num<=pos<=ref_pos+self.max_gap_num:
                            #         check_pos = False
                        if pos-1 == 0 or self.position_to_read[pos] != self.position_to_read[pos-1] or (currK < self.mid_k and pos not in prev_pos):
                        # if pos in self.seq_start or (currK < self.mid_k and pos not in prev_pos):
                            continue
                        if self.pseudogene[pos-1] != 'N':
                            key = self.pseudogene[pos-1]
                            edgeScore[key] += 1
                            if key not in positions:
                                positions[key] = [pos-1]
                            else:
                                positions[key].append(pos-1)
                    for pos in rc_res:
                            # check_pos = True
                            # for ref_pos in rc_prev_pos:
                            #     if ref_pos-self.max_gap_num<=pos+currK-2<=ref_pos+self.max_gap_num:
                            #         check_pos = False
                        if pos + currK-1 >= len(self.position_to_read) or self.position_to_read[pos + currK-2] != self.position_to_read[pos + currK-1] or (currK < self.mid_k and pos+ currK-2 not in prev_pos):
                        # if pos + currK-2 in self.seq_end or (currK < self.mid_k and pos+ currK-2 not in prev_pos):
                            continue
                        if self.pseudogene[pos+currK-1] != 'N':
                            key = complement_base(self.pseudogene[pos+currK-1])
                            edgeScore[key] += 1
                            if key not in positions:
                                positions[key] = [pos+currK-1]
                            else:
                                positions[key].append(pos+currK-1)
                    has_value_exceeding = any(value >= self.mid_threshold for value in edgeScore.values())
                    if currK <= self.mid_k and not has_value_exceeding:
                        # positions = {}
                        if currK > self.min_k:
                            if currK-self.gap_len >= self.min_k:
                                edgeScore = self.suffix_to_base(node,positions,currK-self.gap_len)
                            else:
                                edgeScore = self.suffix_to_base(node,positions,self.min_k)
                        currK = self.min_k
                    # score_list = paths_score_list[node]
                    for base, count in edgeScore.items():
                        if (count >= self.mid_threshold and currK >= self.mid_k) or count >= self.min_threshold:
                            if node in paths_score_list and max_confidence_path != node:
                                del paths_score_list[node]
                            paths_score_list[base+node] = score_list + [(currK/self.kmer_len)*count]
                            if currK < self.kmer_len:
                                tmpPaths[base+node] = [pathsA[node][0]+(currK/self.kmer_len)*count,currK+1]
                                tmpPaths[base+node] = tmpPaths[base+node]+positions[base]
                                # score_list = paths_score_list[node]
                            else:
                                kmer = base + node[:self.kmer_len-1]
                                res,rc_res = self.find_ug_positions(kmer)
                                    # if len(res) + len(rc_res) == 0 or count < self.min_threshold:
                                if len(res) + len(rc_res) == 0:
                                    tmpPaths[base + node] = [pathsA[node][0]+(currK/self.kmer_len)*count,currK]
                                    tmpPaths[base+node] = tmpPaths[base+node]+positions[base]
                                    # score_list = paths_score_list[node]
                                else:
                                    flag = True
                                    # if kmer[:self.kmer_len-1] not in self.unitigs_prefix and DNA_rev_complement(kmer[:self.kmer_len-1]) not in self.unitigs_suffix:
                                    if len(res) != 0:
                                        break_info[base] = (0,res[0])
                                    else:
                                        break_info[base] = (1,rc_res[0])
                                    confidence = self.compute_path_confidence(node,pathsA[node][0],unitig)
                                    if max_confidence == 0:
                                        max_confidence_path = node
                                        max_count = pathsA[node][0]
                                        max_confidence = confidence
                                        max_break_info = break_info
                                    else:
                                        if confidence > max_confidence:
                                            max_confidence_path = node
                                            max_count = pathsA[node][0]
                                            max_confidence = confidence
                                            max_break_info = break_info
                    temp2 = time.time()
                    # self.t3 += temp2 - temp1
                    if len(tmpPaths) > 0 or len(break_info) > 0:
                        break
                    else:
                        currK -= 1
                # if len(break_info) > 0:
                #     confidence = self.compute_path_confidence(node,pathsA[node][0],unitig)
                #     if max_confidence == 0:
                #         max_confidence_path = node
                #         max_count = pathsA[node][0]
                #         max_confidence = confidence
                #         max_break_info = break_info
                #     else:
                #         if confidence > max_confidence:
                #             max_confidence_path = node
                #             max_count = pathsA[node][0]
                #             max_confidence = confidence
                #             max_break_info = break_info
                if currK < self.min_k:
                # if currK < self.mid_k:
                    break_paths[node] = (pathsA[node][0],self.kmer_len)
            if len(tmpPaths) > 0:
                pathsA = tmpPaths
            else:
                best_path = ''
                count = 0
                if len(max_confidence_path) > 0:
                    best_path = max_confidence_path
                    count = max_count
                    for base in max_break_info.keys():
                        kmer =  base + best_path[0:self.kmer_len-1]
                        (fb,kmer_position) = max_break_info[base]
                        self.l_record_break_point(fb,kmer_position)
                    # self.unitigs_prefix.add(best_path[:self.kmer_len-1])
                    # self.rc_unitigs_suffix.append(DNA_rev_complement(best_path[:self.kmer_len-1]))
                else:
                    best_path,count = self.choose_best_path(break_paths,unitig)
                self.unitigs_suffix.discard(unitig[:self.kmer_len-1])
                self.unitigs_prefix.add(best_path[:self.kmer_len-1])
                self.unitigs_ex_info[unitig].left_extension = best_path[:len(best_path)-len(unitig)]
                self.unitigs_ex_info[unitig].left_score = count
                # self.unitigs_ex_info[unitig].left_extension_list = paths_score_list[best_path]
                return best_path[:self.kmer_len-1]
            i += 1
        break_paths.update(pathsA)
        best_path,count = self.choose_best_path(break_paths,unitig)
        self.unitigs_ex_info[unitig].left_extension = best_path[:len(best_path)-len(unitig)]
        self.unitigs_ex_info[unitig].left_score = count
        # self.unitigs_ex_info[unitig].left_extension_list = paths_score_list[best_path]
        return best_path[:self.kmer_len-1]

        
        
        

    def get_max(self,node,count,unitig,max_confidence_break_path,max_break_confidence):
        confidence = self.compute_path_confidence(node,count,unitig)
        if max_break_confidence == 0:
            max_confidence_break_path = node
            max_break_confidence = confidence
        else:
            if confidence > max_break_confidence:
                max_confidence_break_path = node
                max_break_confidence = confidence
        return max_confidence_break_path,max_break_confidence

    
    def choose_best_path(self,break_paths,unitig):
        max_confidence = 0
        best_path = ''
        count = 0
        for path in break_paths.keys():
            if len(path)-len(unitig) == 0:
                return path,break_paths[path][0]
            else:
                confidence = break_paths[path][0]/(len(path)-len(unitig))
            if confidence > max_confidence:
                max_confidence = confidence
                count = break_paths[path][0]
                best_path = path
        return best_path,count

    def compute_path_confidence(self,path,count,unitig):
        if len(path)-len(unitig) == 0:
            return 0
        return count/(len(path)-len(unitig))

    # def check_tail(self, tmpPaths, unitig, flag):
    #     last_kmer_map = {}
    #     first_kmer_map = {}
    #     for path, (count,currK) in tmpPaths.items():
    #         last_kmer = path[len(path)-self.kmer_len:]
    #         first_kmer = path[0:self.kmer_len]
    #         confidence = count/(len(path)-len(unitig))  #计算路径置信度
    #         if last_kmer not in last_kmer_map or confidence > last_kmer_map[last_kmer][1]:
    #             last_kmer_map[last_kmer] = (path, (count,currK))
    #         if first_kmer not in first_kmer_map or confidence > first_kmer_map[first_kmer][1]:
    #             first_kmer_map[first_kmer] = (path, (count,currK))
    #     tmpPaths.clear()  # 清空 tmpPaths
    #     if flag == 0:
    #         for path, (count,currK) in last_kmer_map.values():
    #             tmpPaths[path] = (count,currK)
    #     else:
    #         for path, count in first_kmer_map.values():
    #             tmpPaths[path] = (count,currK)
    #     return tmpPaths
        
    

    def change_order_right1(self,currK):
        tmpNodes = []
        link_unitigs = []
        for node in self.currNodes:
            edgeScore = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
            currK = self.kmer_len if len(node) > self.kmer_len else len(node)
            pNode = node[len(node)-currK:]
            temp1 = time.time()
            result = self.fmi.get_range(pNode[1:],0)
            res = self.fmi.get_offsets(result[0],result[1])
            temp2 = time.time()
            self.findright2 += temp2-temp1
            if res is None:
                break
            if len(res)-1>self.copy and currK < 11:
                self.findleft1 += 1
            seq_ids={}
            poss={}
            for pos in res:
                seq_id = self.position_to_sequence.get(pos,None)
                start,end = self.sequence_info.get(seq_id)
                if start <= pos < end and pos + len(pNode) - 1 < end-1:
                    # edgeScore[self.pseudogene[pos+res[0]:pos+res[0]+1]] += 1
                    key = self.pseudogene[pos + len(pNode) - 1:pos + len(pNode)]
                    # 更新 edgeScore 计数
                    edgeScore[key] += 1
                    if key not in seq_ids:
                        seq_ids[key] = []
                    if key not in poss:
                        poss[key] = []
                    seq_ids[key].append(seq_id)
                    poss[key].append(pos + len(pNode))
            for block in edgeScore.keys():
                if 2 <= edgeScore[block]:
                    if currK < self.kmer_len:
                        tmpNode = self.extend_node(seq_ids[block],poss[block],node[len(node)-currK:len(node)]+block,2,False)
                    else:
                        tmpNode[node[len(node)-currK+1:len(node)]+block]=edgeScore[block]
                    for temp in tmpNode.keys():
                        link_unitigs = self.check_link(node,temp,1)
                        if len(link_unitigs)==0:
                            tmpNodes.append(temp)
                            self.graph.add_node(temp)
                            if currK < self.kmer_len:
                                self.graph.add_edge(node,temp,coverage=tmpNode[temp],overlap=len(pNode))
                            else:
                                self.graph.add_edge(node,temp,coverage=tmpNode[temp],overlap=self.kmer_len-1)
                        else:
                            for link_unitig in link_unitigs:
                                self.graph.add_edge(node,link_unitig,overlap=currK-1)
        return tmpNodes

    # def find_right_paths(self,currK):
    #     tmpNodes = []
    #     # temp1 = time.time()
    #     if(currK == self.kmer_len):
    #         tmpNodess = self.change_order_right(currK)
    #         if(len(tmpNodes)==0):
    #             currK -= 1
    #         else:
    #             return tmpNodes,currK
    #     if(self.min_len <= currK < self.kmer_len):
    #         while(len(tmpNodes) == 0 and currK >= self.min_len-1):
    #             tmpNodes = self.change_order_right(currK)
    #             if len(tmpNodes) > 0 or len(link_unitigs) > 0:
    #                 currK += 1
    #                 return tmpNodes,currK
    #             currK -= 1
    #     return tmpNodes,currK

    def score_bases(self, curr_bases, indices, positions, threshold, indices_idx, isleft):
        bases_above_threshold = {}
        for idx in indices_idx:
            read = self.reads[indices[idx]]
            pos = positions[idx]
            if not isleft:
                positions[idx] = positions[idx] + 1  # 将位置后移
            else:
                positions[idx] = positions[idx] - 1
            if 0 <= pos < len(read):
                base = read[pos]
                # 无论是否超过阈值，都记录 idx
                if base not in bases_above_threshold:
                    bases_above_threshold[base] = []  # 初始化一个列表来记录 idx
                bases_above_threshold[base].append(idx)
                
                # 增加碱基的计数
                if base in curr_bases:
                    curr_bases[base] += 1
        
        # 过滤出超过阈值的碱基，保留所有记录的 idx
        bases_above_threshold = {base: idx_list for base, idx_list in bases_above_threshold.items() if curr_bases[base] >= threshold}
        return bases_above_threshold


    def score_next_bases(self):
        self.nextBases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        maxScore = 0
        for m in self.nextBases.keys():
            self.nextBases[m] = self.score_next_base(m)
            if self.nextBases[m] > maxScore:
                maxScore = self.nextBases[m]
        return maxScore

    def score_prev_bases(self):
        self.prevBases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        maxScore = 0
        for m in self.prevBases.keys():
            self.prevBases[m] = self.score_prev_base(m)
            if self.prevBases[m] > maxScore:
                maxScore = self.prevBases[m]
        return maxScore

    def score_next_base(self, aBase):
        pNode = self.pNode['a']
        nKmer1 = pNode[1:len(pNode)] + aBase[0:1]
        if not self.kmers.get(nKmer1):
            return 0
        return self.kmers.get(nKmer1)

    def score_prev_base(self, aBase):
        pNode = self.pNode['b']
        nKmer1 = aBase[0:1] + pNode[0:len(pNode) - 1]
        if not self.kmers.get(nKmer1):
            return 0
        return self.kmers.get(nKmer1)
    
    def remove_low_cov_kmers(self, min_cov=10):
        if min_cov > 0:
            nKmlist = {}
            for a in self.kmers:
                if (self.kmers[a] > min_cov):
                    nKmlist[a] = self.kmers[a]
            self.kmers = nKmlist
