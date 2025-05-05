import random
import math
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from random import choices
from typing import Dict, List, Tuple, Set

def generate_paired_reads(fasta_file, fastq_file_1, fastq_file_2, depth, fragment_size=300, read_length=150):

    def simulate_quality(length):
        """生成模拟的高质量分数。"""
        return ''.join(['~'] * length)  # 使用最高质量 '~' 代表质量分数

    with open(fastq_file_1, 'w') as fq_out1, open(fastq_file_2, 'w') as fq_out2:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            seq_length = len(sequence)

            # 计算生成的片段数量
            num_fragments = (depth * seq_length) // fragment_size

            for i in range(num_fragments):
                # 随机选择一个起始位置
                start = random.randint(0, seq_length - fragment_size)
                fragment = sequence[start:start + fragment_size]

                # 在片段中添加随机错误
                # e_DNA = [fragment]
                # e_DNA = add_rand_sub_new(e_DNA, 0.005)  # 替换错误
                # e_DNA = add_rand_del_new(e_DNA, 0.0025)  # 删除错误
                # e_DNA = add_rand_indel_new(e_DNA, 0.0025)  # 插入错误

                # # 获取带错误的片段
                # fragment = e_DNA[0]

                # 生成配对端的两个读段
                read1 = fragment[:read_length]
                read2 = Seq(fragment[-read_length:]).reverse_complement()

                # 模拟读段的质量分数
                quality1 = simulate_quality(read_length)
                quality2 = simulate_quality(read_length)

                # 写入 FASTQ 文件
                fq_out1.write(f"@{record.id}_{i}/1\n{read1}\n+\n{quality1}\n")
                fq_out2.write(f"@{record.id}_{i}/2\n{read2}\n+\n{quality2}\n")

def concatenate_fasta_to_single_line(fasta_file):
    """
    将FASTA文件中的所有序列拼接成一行，并覆盖写回原FASTA文件。
    """
    concatenated_sequence = []

    # 读取FASTA文件并拼接序列
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()  # 移除换行符和空白符
            if not line.startswith('>'):  # 忽略描述行
                concatenated_sequence.append(line)  # 收集序列部分

    # 将拼接后的序列合并为一行
    sequence_in_one_line = ''.join(concatenated_sequence)

    # 重新写入原文件，保留描述行，并将序列写为一行
    with open(fasta_file, 'w') as f:
        f.write(">concatenated_sequence\n")  # 写入新的描述行
        f.write(sequence_in_one_line + '\n')  # 写入拼接后的序列

def	charN(str, N):
	if N < len(str):
		return str[N]
	return 'X'


def oxor(str1, str2):
	length = max(len(str1), len(str2))
	return ''.join(chr(ord(charN(str1, i)) ^ ord(charN(str2, i))) for i in range(length))

def	xor(str1, str2):
	length = max(len(str1), len(str2))
	if len(str1) > len(str2):
		str2 = str2 + bytes(len(str1) - len(str2))
	if len(str2) > len(str1):
		str1 = str1 + bytes(len(str2) - len(str1))
	allBytes = b''
	i = 0
	while i < length:
		allBytes = allBytes + bytes([str1[i] ^ str2[i]])
		i =  i + 1
	return allBytes


def randChunkNums(degree, num_chunks):
	size = random.randint(1,max(5, int(num_chunks*3/5)))
	return random.sample(range(num_chunks), size)


# def get_degrees(N, k, symbol_index=-1, distribution_name='robust'):
# 	print(distribution_name)
# 	if distribution_name == "ideal":
# 		probabilities = ideal_distribution(N)
# 	elif distribution_name == "robust":
# 		probabilities = robust_distribution(N, 0.01, 0.01)
# 	else:
# 		probabilities = None
#
# 	population = list(range(0, N))
# 	if symbol_index > 0:
# 		random.seed(symbol_index)
# 	#Not initialized yet
# 	else:
# 		return -1
# 	return choices(population, probabilities, k=k)



def get_degrees(N, k, symbol_index=-1, delta=0.01, c_value=0.01, distribution_name='robust'):
#def get_degrees(N, k, symbol_index=-1, delta=0.1, c_value=0.2, distribution_name='robust'):
	#print(distribution_name)
	if distribution_name == "ideal":
		probabilities = ideal_distribution(N)
	elif distribution_name == "robust":
		probabilities = robust_distribution(N, c_value, delta)
	else:
		probabilities = None

	population = list(range(0, N))
	if symbol_index > 0:
		random.seed(symbol_index)
	#Not initialized yet
	else:
		return -1
	return choices(population, probabilities, k=k)


def generate_chunk_nums(blocks_quantity, degree, symbol_index):
	random.seed(symbol_index)
	indexes = random.sample(range(blocks_quantity), degree[0])
	return indexes


def ideal_distribution(K):

	probabilities = [0, 10 / K]
	probabilities += [1 / (k * (k - 1)) for k in range(2, K)]
	probabilities_sum = np.sum(probabilities)
	probabilities /= probabilities_sum

	# assert probabilities_sum >= 1 - epsilon and probabilities_sum <= 1 + epsilon, "The ideal distribution should be standardized"
	return probabilities

def robust_distribution(K, c_value = 0.5, robust_failure_probability=0.1):

	# print('c value=', end='\t')
	# print(c_value, end='\t')
	# print('epsilon', end = '\t')
	# print(robust_failure_probability)

	S = c_value * math.log(K / robust_failure_probability) * math.sqrt(K)
	# print(S)
	# M = round(K/S)
	M = round(K / S)
	# print(M)
	extra_proba = [0] + [1/(i * M) for i in range(1, M-1)]
	extra_proba += [S * math.log(S / robust_failure_probability) / K]  # Spike at M
	extra_proba += [0 for k in range(M, K)]

	probabilities = np.add(extra_proba, ideal_distribution(K))
	# print(np.sum(probabilities))
	probabilities /= np.sum(probabilities)
	#probabilities_sum = np.sum(probabilities)
	#assert probabilities_sum >= 1 - epsilonc and probabilities_sum <= 1 + epsilonc, "The robust distribution should be standardized"
	return probabilities

def bytesToDNA(manyBytes):
	dnastr = ''
	for aByte in manyBytes:
		dnastr = dnastr + byteToDNA(aByte)
	return dnastr

def DNAToBytes(dnaStr):
	i = 0
	nBytes=b''
	while i < len(dnaStr):
		nBytes = nBytes + simDNAToByte(dnaStr[i:i+4])
		i = i + 4
	return nBytes

def byteToDNA(abyte):
	#define the converter of four bits to DNA chars
	converter = ('AT', 'AG', 'AC', 'AA',  'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC',   'CC', 'CT', 'CA', 'CG')
	octNum = abyte
	#octNum = ord(abyte)
	dnastr = converter[octNum % 16]
	octNum = octNum//16
	dnastr = converter[octNum % 16] + dnastr
	return dnastr



def simDNAToByte(dnaStr):
	converter = ('AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC', 'CC', 'CT', 'CA', 'CG')
	bytesConverter = {}
	i = 0
	while i < 16:
		j = 0
		while j < 16:
			fourBases = converter[i] + converter[j]
			bytesConverter[fourBases] = bytes([i*16 + j])
			j = j + 1
		i = i + 1
	return bytesConverter[dnaStr[0:2] + dnaStr[2:4]]

def DNAToByte(dnaStr):
	converter = ('AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC', 'CC', 'CT', 'CA', 'CG')
	bytesConverter = {}
	i = 0
	while i < 16:
		j = 0
		while j < 16:
			fourBases = converter[i] + converter[j]
			bytesConverter[fourBases] = bytes([i*16 + j])
			j = j + 1
		i = i + 1
	return bytesConverter[dnaStr[0:2] + dnaStr[3:5]]




def compressDNA(dnaStr):
	#return dnaStr
	converter = ('A', 'T', 'G', 'C')
	bytesConverter = {}
	i = 0
	while i < 4:
		j = 0
		while j < 4:
			k = 0
			while k < 4:
				m = 0
				while m < 4:
					fourBases = converter[i] + converter[j] + converter[k] + converter[m]
					bytesConverter[fourBases] = bytes([i*4*4*4 + j*4*4 + k*4 + m])
					m = m + 1
				k = k + 1
			j = j + 1
		i = i + 1
	i = 0
	nBytes = b''
	while i < len(dnaStr):
		nBytes = nBytes + bytesConverter[dnaStr[i:i+4]]
		i = i + 4
	return dnaStr

def depressDNA(manyBytes):
	#define the converter of four bits to DNA chars
	converter = ('A', 'T', 'G', 'C')
	dnastr = ''
	for aByte in manyBytes:
		fourBase = ''
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		dnastr = dnastr + fourBase
	return dnastr



#2019-5-14
def calc_gc(dnastr):
	dnastr = dnastr.upper()
	gc_num = dnastr.count("G") + dnastr.count("C")
	return gc_num/len(dnastr)

def max_homo_len(dnastr):
	dnastr = dnastr.upper()
	max_len = 0
	pre_len = 0
	last_c = ''
	for c in dnastr:
		if c == last_c:
			pre_len = pre_len + 1
		else:
			if pre_len > max_len:
				max_len = pre_len
			pre_len = 1
		last_c = c
	if pre_len > max_len:
		max_len = pre_len
	return max_len

#2019-05-15
#2020-03-08
#def check_dna(dnastr, min_gc=0.45, max_gc=0.55, max_homo=5):
def check_dna(dnastr, min_gc=0.45, max_gc=0.55, max_homo=5):
	# print(dnastr)
	gc_rate = calc_gc(dnastr)
	if gc_rate > max_gc:
		return False
	if gc_rate < min_gc:
		return False
	homo_poly_len = max_homo_len(dnastr)
	if homo_poly_len > max_homo:
		return False
	return True

#2019-05-16
def rev_seq(dnastr):
	complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
	dnastr = dnastr.upper()
	rev_dna = ""
	for i in dnastr:
		rev_dna += complement[i]
	rev_dna = rev_dna[::-1]
	return rev_dna

#20190625
def num_randint(a,b,num):
	i = 0
	int_nums = []
	while i < num:
		int_nums.append(random.randint(a, b))
		i += 1
	return int_nums


def randomATGC():
	a = "ATGC"
	ar = random.randint(0, 3)
	return a[ar:ar + 1]

def randomATGC(noATGC='', len=1):
        all_bases = ["A", "G", "C", "T"]
        allow_ATGC = []

        for aBase in all_bases:
            if not aBase in noATGC:
                allow_ATGC.append(aBase)

        random_sequence = ''
        for i in range(0, len):
            random_sequence += random.choice(allow_ATGC)
        return random_sequence


def randomDNA(len):
	i = 1
	dna = ""
	while(i <=len):
		dna= dna + randomATGC()
		i = i + 1

	return dna

def dna_to_numbers(sequence):
    # 创建一个字典来映射每个碱基到其对应的数字
        base_to_number = {'A': 1, 'G': 2, 'C': 3, 'T': 4}

        # 使用列表生成式来生成对应的数字数组
        number_sequence = [base_to_number[base] for base in sequence]

        return number_sequence

def complement_base(base):
    # 定义碱基的互补关系
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return complement.get(base.upper())

def DNA_complement(sequence):
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	return sequence.upper()

def DNA_reverse(sequence):
	sequence = sequence.upper()
	return sequence[::-1]

def DNA_rev_complement(sequence):
	sequence = DNA_complement(sequence)
	return DNA_reverse(sequence)

def expntl(L):
	"""
    negative exponential distribution
    return a double random number, L is the mean value
    """
	u = random.randint(0,2**L)
	return int(np.math.log(u)/math.log(2))



def file_to_array(file):
	f = open(file)
	arr = []
	line = f.readline()
	while line.strip():
		arr.append(line.strip())
		line = f.readline()
	return arr


def kmers_of_str(str, kmer_len=21,step_len=1):
	kmers = {}
	i = 0
	if len(str) >= kmer_len:
		i = 0
		kmstr = ''
		while i <= len(str) - kmer_len:
			kmstr = str[i:i + kmer_len]
			kmers[kmstr] = 1
			i = i + step_len
	return kmers.keys()

def seq_to_fastq(sequences, output_file):
    with open(output_file, 'w') as f:
        for index, sequence in enumerate(sequences):
            seq_id = f"seq{index}"  # 使用序列的下标作为 id
            quality_scores = 'I' * len(sequence)  # 最高质量分数

            # 写入 FASTQ 格式内容
            f.write(f"@{seq_id}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality_scores}\n")

def calculate_length_sum(file_path):
    total_length = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                # 找到包含length的部分
                length_part = line.split('length=')  # 以 'length=' 分割
                if len(length_part) > 1:
                    length_value = length_part[1].strip()  # 获取长度值并去掉空白字符
                    total_length += int(length_value)  # 累加长度值
    return total_length

def read_fastq_sequences(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        while True:
            # 读取四行
            header = file.readline().strip()  # 第一行：序列的描述信息
            if not header:  # 文件结束
                break
            sequence = file.readline().strip()  # 第二行：实际的序列
            file.readline().strip()  # 第三行：'+'
            file.readline().strip()  # 第四行：质量值
            
            # 将序列添加到列表中
            sequences.append(sequence)
    return sequences

def get_contig(d: Dict[str, int], x: str) -> Tuple[str, int]:
    """从字典中获取 contig 和它的覆盖度"""
    contig = x
    coverage = d[x]
    return contig, coverage

def fw(kmer: str) -> List[str]:
    """生成k-mer的前缀，用于查找匹配"""
    return [kmer[1:] + base for base in "ATCG"]

def all_contigs(d: Dict[str, int], k: int) -> Tuple[Dict[int, Tuple[List[Tuple[int, str]], List[Tuple[int, str]]]], List[str]]:
    """
    构建 contig 之间的连接关系，并返回 contig 序列和关系图。
    G[i][0] 代表 contig i 的尾部连接关系，G[i][1] 代表头部连接关系。
    """
    done: Set[str] = set()
    r: List[str] = []  # 保存所有的 contig 序列

    # 遍历字典中的 contig 并记录到集合中
    for x in d:
        if x not in done:
            s, _ = get_contig(d, x)
            for y in [x, DNA_rev_complement(x)]:
                done.add(y)
            r.append(s)

    G = {}  # 用于存储 contig 的连接关系
    heads = {}  # 保存 contig 的头部 k-mer 到 contig 的映射
    tails = {}  # 保存 contig 尾部 k-mer 的反向互补映射

    # 初始化 G，并将头部和尾部 k-mer 记录在 heads 和 tails 中
    for i, x in enumerate(r):
        G[i] = ([], [])  # 每个 contig 的头部和尾部的连接列表
        heads[x[:k]] = (i, '+')  # 头部的 k-mer
        tails[DNA_rev_complement(x[-k:])] = (i, '-')  # 尾部的反向互补 k-mer

    # 建立 contig 之间的连接关系
    for i, x in enumerate(r):
        # 处理 contig 尾部的 k-mer 匹配关系
        for y in fw(x[-k:]):
            if y in heads:
                G[i][0].append(heads[y])  # 尾部连接到另一个 contig 的头部
            if y in tails:
                G[i][0].append(tails[y])  # 尾部连接到另一个 contig 的尾部

        # 处理 contig 头部的 k-mer 匹配关系
        for z in fw(DNA_rev_complement(x[:k])):
            if z in heads:
                G[i][1].append(heads[z])  # 头部连接到另一个 contig 的头部
            if z in tails:
                G[i][1].append(tails[z])  # 头部连接到另一个 contig 的尾部

    return G, r

def print_fasta_format_with_links_to_file(G, contigs,k, d, filename):
    """将 contigs 和连接关系按指定格式输出到文件"""
    with open(filename, 'w') as f:
        for i, contig in enumerate(contigs):
            # 输出 contig 的信息，包括长度、覆盖度和 k-mer 覆盖率
            coverage = float(f"{d[contig]:.2f}")  # 保留两位小数的覆盖率
            count = round(coverage * (len(contig)-k+1))  # 计算 KC

            # 构建 header
            header = f">contig_{i} LN:i:{len(contig)} KC:i:{count} km:f:{coverage}"

            # 添加 contig 的尾部连接关系
            for j, direction in G[i][0]:  # 尾部连接
                header += f" L:+:{j}:{direction}"

            # 添加 contig 的头部连接关系
            for j, direction in G[i][1]:  # 头部连接
                header += f" L:-:{j}:{direction}"

            # 写入文件，添加换行符
            f.write(header + "\n")
            f.write(contig + "\n")


def seq_to_fasta(sequences, output_file, prefix="seq"):
    # 打开输出文件进行写入
    with open(output_file, 'w') as f:
        for i, sequence in enumerate(sequences):
            # 如果序列长度不足1024，则补充N到1024长度
            if len(sequence) < 1024:
                sequence = sequence + 'N' * (1024 - len(sequence))
            # 自动生成标识符
            identifier = f"{prefix}{i+1}"
            # 写入FASTA格式
            f.write(f">{identifier}\n")
            f.write(f"{sequence}\n")




#2020-05-03
def kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer not in dict:
			return False
	return True

def any_kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer in dict:
			return True
	return False

def read_file(file):
	file1 = open(file, 'rb')
	filebytes = file1.read()
	file1.close()
	return filebytes




# def fasta_to_dict(fasta_file):
#     unitigs = []

#     # 读取 FASTA 文件
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         unitig = str(record.seq)
#         unitigs.append(unitig)

#     return unitigs

def fasta_to_dict(fasta_file):
    unitigs = {}

    # 读取 FASTA 文件
    for record in SeqIO.parse(fasta_file, "fasta"):
        unitig = str(record.seq)
        
        # 提取 `km:f:` 后的值
        description = record.description
        km_value = None
        if "km:f:" in description:
            try:
                # 获取 "km:f:" 后面的值
                km_value = float(description.split("km:f:")[1].split()[0])
            except ValueError:
                print(f"Error parsing km value in record: {description}")
        
        if km_value is not None:
            unitigs[unitig] = km_value

    return unitigs

def read_sequences_to_array(fasta_file):
    """
    读取FASTA文件，将所有序列存入数组。
    :param fasta_file: str, FASTA文件路径
    :return: list, 包含所有序列的数组
    """
    sequences = []  # 用于存储所有序列
    current_sequence = ""  # 存储当前序列
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # 如果是标题行，跳过
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""  # 重置当前序列
            else:
                # 拼接序列行
                current_sequence += line

        # 将最后一个序列存入数组
        if current_sequence:
            sequences.append(current_sequence)

    return sequences

def read_fasta_to_list(fasta_file):
    """
    读取 FASTA 文件并将每个序列存储到一个列表中。

    :param fasta_file: FASTA 文件的路径。
    :return: 包含所有序列的列表。
    """
    sequences = []
    current_sequence = []

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(''.join(current_sequence))
                    current_sequence = []
            else:
                current_sequence.append(line)
        if current_sequence:
            sequences.append(''.join(current_sequence))

    return sequences
