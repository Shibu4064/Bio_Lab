from Bio import Entrez,SeqIO
Entrez.email="hrithikworks22gmail.com"

handle=Entrez.einfo(db='pubmed')
result=Entrez.read(handle)
handle.close()
#result
print(result.keys())

handle=Entrez.efetch(db='nucleotide',id='NM_001126.2',rettype='fasta',retmode='text')
record=SeqIO.read(handle,'fasta')
print(f'''Sequence ID:{record.id},Sequence Name:{record.name},Sequence Length:{len(record)},
Sequence Description:{record.description}''')

SeqIO.write(record,'testseq.fasta','fasta')

# Fasta file reading code

record=SeqIO.read('/content/testseq.fasta','fasta')
print(f'''Sequence ID:{record.id},Sequence Name:{record.name},Sequence Length:{len(record)},
Sequence Description:{record.description}''')

handle=Entrez.efetch(db='nucleotide',id=['34577062','186972394'],rettype='fasta',retmode='text')
records=list(SeqIO.parse(handle,'fasta'))
print(len(records))
print(records[0].id,len(records[0]))
print(records[1].id,len(records[1]))

SeqIO.write(records,'MultiSeq.fasta','fasta')

records=SeqIO.parse('/content/MultiSeq.fasta','fasta')
for record in records:
  print(record.id,len(record))


# 1A: Compute the Number of Times a Pattern Appears in a Text

from Bio.Seq import Seq
seq=Seq('GCGCG')
seq.count_overlap('GCG')

# 1D: Find all Occurences of a Pattern in a String

seq='GATATATGCATATACTT'
kmer='ATAT'
k=len(kmer)
i=0
while(i+k<len(seq)):
    if(seq[i:i+k]==kmer):
        print(i,end=" ")
    i=i+1

# 1F: Find a Position in a Genome Minimizing the Skew

import numpy as np
import matplotlib.pyplot as plt
def skew(dna):
    map={'A':0,'T':0,'G':1,'C':-1}
    value=np.array([0])
    value=np.append(value,[map[x] for x in dna],axis=0)
    y=np.cumsum(value)
    plt.figure(figsize=(15,10))
    plt.plot(y)
    plt.xticks(range(len(dna)),dna)
    plt.grid()
    plt.show()
    return np.where(y==y.min())
dna='CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
res=skew(dna)
print(" ".join(map(str,res)))

# 1J: Find Frequent Words with Mismatches and Reverse Compliments

seq='ACGTTGCATGTCGCATGATGCATGAGAGCT'
k=4
d=1
dic={}
mx=0
i=0
complements=dict(zip('ATGC','TACG'))
def getrc(s):
    return ''.join([complements[base] for base in s][::-1])
def create(s):
    if len(s)==k:
        dic[s]=0
        return
    create(s+'A')
    create(s+'T')
    create(s+'G')
    create(s+'C')
create('')
#print(dic)
for kmer in dic:
    i=0
    while i+k<=len(seq):
        hd=0
        for j in range(k):
            if seq[i+j]!=kmer[j]:
                hd=hd+1
        if hd<=d:
            dic[kmer]=dic[kmer]+1
        i=i+1
for kmer in dic:
    mx=max(mx,dic[kmer]+dic[getrc(kmer)])
for kmer in dic:
    if dic[kmer]+dic[getrc(kmer)]==mx:
        print(kmer,end=" ")

# 1K: Generate the Frequency Array of a String

def generate_frequency_array(s,k):
    frequency_array=[0]*(4**k)
    base_to_num={'A':0,'C':1,'G':2,'T':3}
    for i in range(len(s)-k+1):
        kmer=s[i:i+k]
        num=0
        for j in range(k):
            num=num*4+base_to_num[kmer[j]]
        frequency_array[num]+=1
    return frequency_array
dna_string="ACGCGGCTCTGAAA"
k=2
frequency_array=generate_frequency_array(dna_string,k)
print(" ".join(map(str,frequency_array)))

plt.bar(range(len(frequency_array)),frequency_array)
plt.xlabel('K mer index')
plt.ylabel('frequency')
plt.show()

# 1L: Implement Pattern to Number

dic={'A':0,'C':1,'G':2,'T':3}
seq='AGTG'
i=len(seq)-1
s=0
f=1
while i>=0:
    s=s+f*dic[seq[i]]
    f=f*4
    i=i-1
print(s)

# 1M: Implement Number To Pattern

dic=dict(zip(range(4),'ACGT'))
n=45
k=4
i=k
f=pow(4,k-1)
while k>0:
    print(dic[int(n/f)],end='')
    n=n%f
    f=f/4
    k=k-1

# 1N: Generate The d-Neighborhood of a String

import itertools
def hamming_distance(s1,s2):
    return len(s1)-sum([s1[i]==s2[i] for i in range(len(s1))])
def dNeighborhood(dna,d):
    kmers=[''.join(kmer) for kmer in itertools.product('ACGT',repeat=len(dna))]
    return [kmer for kmer in kmers if hamming_distance(dna,kmer)<=d]
dna='ACG'
d=1
res=dNeighborhood(dna,d)
print(" ".join(map(str,res)))

# 2B: Find a Median String

from itertools import product
def hamming_distance(s1, s2):
    return sum(1 for x, y in zip(s1, s2) if x != y)
def find_median_string(k, sequences):
    min_distance = float("inf")
    median = ""
    for pattern in product("ACGT", repeat=k):
        pattern = ''.join(pattern)
        total_distance = sum(min(hamming_distance(pattern, sequence[i:i+k]) for i in range(len(sequence) - k + 1)) for sequence in sequences)
        if total_distance < min_distance:
            min_distance = total_distance
            median = pattern
    return median
def read_input(filename):
    with open(filename, 'r') as file:
        k = int(file.readline().strip())
        sequences = [line.strip() for line in file]
    return k, sequences
input_filename = '/content/rosalind_ba2b.txt'
k, sequences = read_input(input_filename)
median_string = find_median_string(k, sequences)
print(median_string)

# 2C: Find a Profile Most Probable K mer in a String

f=open("/content/rosalind_ba2c.txt")
seq=f.readline()
seq=seq.replace('\n','')
k=int(f.readline())
profile=[]
for lis in f:
  profile.append([float(i) for i in lis.split(" ")])
#print(profile)
kmers=[]
i=0
dic={}
ids=dict(zip("ACGT",range(4)))
while i+k<=len(seq):
  dic[seq[i:i+k]]=0
  i=i+1
res=seq[0:k]
for kmer in dic:
  pro=float(1)
  for i in range(k):
    #print(kmer)
    pro*=profile[ids[kmer[i]]][i]
  dic[kmer]=pro
for kmer in dic:
  if dic[kmer]>dic[res]:
    res=kmer
print(res)

# 2F: Implement Randomized Motif Search

import math
import random
def hamming_distance(text1, text2):
  cnt = 0
  for i in range(len(text1)):
    if not(text1[i] == text2[i]):
      cnt += 1
  return cnt
def score(motifs, consensus):
  total_score = 0
  for motif in motifs:
    total_score += hamming_distance(motif, consensus)
  return total_score
def consensus_string(motifs):
  k = len(motifs[0])
  c_s = ""
  for i in range(k):
    dic = {
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0
    }
    for motif in motifs:
      ch = motif[i]
      dic[ch] += 1
    mx = max(dic.values())
    for key in dic:
      if mx == dic[key]:
        c_s += key
        break
  return c_s
def get_profile_most_probable_kmer(profile, dna, k, t):
  frequency_matrix = profile.copy()
  for key in frequency_matrix:
    frequency_matrix[key] = [(val+1)/(t+4) for val in frequency_matrix[key]]
  n = len(dna)
  kmer_list = []
  max_score = -1.0
  for i in range(n-k+1):
    score = 1.0
    kmer = dna[i:i+k]
    for j, ch in enumerate(kmer):
      score *= frequency_matrix[ch][j]
    max_score = max(max_score, score)
    kmer_list.append((kmer, score))
  for kmer in kmer_list:
    if max_score == kmer[1]:
      return kmer[0]
def create_profile(motifs, k):
  matrix = {
      'A': [0]*k,
      'C': [0]*k,
      'G': [0]*k,
      'T': [0]*k
  }
  for motif in motifs:
    for i, ch in enumerate(motif):
      matrix[ch][i] += 1
  return matrix
def randomized_motif_search_with_pseudocount(dnas, k, t):
  best_score = math.inf
  n = len(dnas[0])
  motifs = []
  for dna in dnas:
    idx = random.randint(0, n-k)
    motifs.append(dna[idx:idx+k])
  best_motifs = motifs.copy()
  while True:
    profile = create_profile(motifs, k)
    motifs.clear()
    for j in range(len(dnas)):
      profile_most_probable_kmer = get_profile_most_probable_kmer(profile, dnas[j], k,  t)
      motifs.append(profile_most_probable_kmer)
    consensus = consensus_string(motifs)
    score_motif = score(motifs, consensus)
    if score_motif < best_score:
      best_score = score_motif
      best_motifs = motifs
    else:
      return best_motifs, best_score
if __name__ == "__main__":
  k, t = 8, 5
  dnas = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
  ]

  # with open('rosalind_ba2f.txt') as file:
  #   f = file.read().strip().split("\n")
  #   first = f[0].split()
  #   k = int(first[0])
  #   t = int(first[1])
  #   dnas = f[1:]
  best_motifs, best_score = randomized_motif_search_with_pseudocount(dnas, k, t)
  for i in range(1,1000):
    motif, _score = randomized_motif_search_with_pseudocount(dnas, k, t)
    if _score < best_score:
      best_motifs = motif
      best_score = _score
  for motif in best_motifs:
    print(motif)

# 2H: Implement Distance Between Pattern and Strings

pattern="CGGCGAT"
texts="TGTGTTGATGGGATGAAGGTTAGGTCGCTGAAGCGTGTACATTAGAGCGTTCTACCCTCGTGGGAAGGTGATACCCAGGAGGGAGTGCAC AGTGAATGTGGCTGTGTGTTTTGTTGAATAAGGGCTCCCACTCGTAGGACATCTAGCACTCGCTCTCGCGGCGCACCGAGTGAGTCTGGT TTGATAGCACCAATCGCTCGGATTGACATTGTCCGTAAGGCTAGCAACATATACGAGGAAACATCGAATGCTTCGCAGTCCTCCATGTAA CCAGCGCGCGAGTATGATTGTCACGGGCCCTCTGCAATGCATCACGGATTTCACATGGCAAGGGCTCAAGAGTCTATGCGACCCCCTTAC ACTATGCCCCCTGTTACGCCTTGTCGCGATAAACAAGCTATTATGGGGTATCTACATTGGACCGTTTGTCCAACAACGGGTCAAGCCCCC AAGCACCGCGCGAAACGCAGTTACTGTGGCCCTAAACACGACTGTCTGCGGGGTATAACTTTTCCTCCTGAGTCCCCTGGAAGAGGCGAA CCCGTGTTCCCCGAAATGGAATACCCGGCGGAGTTCAGTAAACTGTTTTCTCGGTACGTACCAGACACTATTTCCACGCACCCGCCAAAC TGGTATTCTTTAATCGGCAGCTCCAGAGAGCGGCGCTTCGCCCTACCAATAGCAGATATACCACGGTGTCTACGGTTAACCTACCGGGTC GGCCTGCTGTGTCCCGAAAACCGGCGTTGTGACTAAACCATCGATCGTCATCTGTGAGCGTGCTCCATATCATACATCAATAGACTCAGA ATTGCAGGGCTCATCAAAAATAACTAATTGTGACTGAGTGATCTTGGTCCTTGAACGGACCTATAAGTGCGCCCGCGTCCTTCCTCGACC CTAATCATCTCATATTTAAAACTAGTATAACTAGTAAATCCTCAACAGAGCTTTGATGGTCGATTCAGAACCAGCCTTTGGCCACCCACT AGAGTCCAAACTCTCTCAATTGACAAACTTTATTGTATGACCAGATAGACACGGCCCTAGATGGCAGTACTTGCTGTGAGTCAGACGAGA AATGTTGGCGTGACTTCGGAAGAAGACGTCACCGGGAGCATGCTACAAAAATATAAAGCAGCCGAGCGGCCCGGGCAGGTGTCCTGAGGG CCGCACTAAGTTTAGCATGAGCGGTTAATCGCCTGTCCCCAGACACTACTTATCCTAGTATCCATACGATGGACGGAATATCTAAAAGGG CCTTAGCGTCGACTGGGCGCTGCGAATGGCGTTATCTAAGTAACAAGTATCTAGGGTGATGGGAGGACGGGGCTCTCGTTTCGCTCTAGT CTCCATCGGAGAGCATACCTTGTGCCCCTCGTAGTGACCGGTTGGTTTATAGCCGGTGACACCAATCCTCGGAAAAGGAATCGGCGCCCC CATTGGAATTGTGTGCTTAGGTCATTTATCTAGCAGGGTATAGTCTGCACATGAGTCTGATTTGGGAGGCAATGTGTTAGGCCTTACAGG TAATCGTCGCGGGAAGGAAGCGCTGAGCGCAACTAATCGTGCGGTGGGAGCGGACTACGAAGATGGTCTCCTCGGCGCTTTAACCCATAC TAACATCAACGGACTGGTGCTGCCAAGCAGTATGGGCTACGACATATTGAAATACTAGCACTGTTCGCTACCCATGGCGACGGGTCTGCT AAAGGTGCTCCGGTAAACCTATTTCCGACGGGGACGGCAGATCTAATACATAACAGTGACTCTAGTGTTCACGCATTCAACGGACGCCTG TTTCTATTATGGTAAGAATATTGAGCGAACCTGTAACGTGAGAAGCCGTAGTCGCAGTCGTGCCTCGCGAGTGCAGACCGAGAGCCAAGC TCATAGGCTCGGGAGAGATTGACACGAGTGTTCCGGTTATTTTAAGACCGAGAGAGCAGTAGACATAGTGCCTGTGCTTCGTTGCACATA AACGGGTCACCAGAAGTAGGAAATCCCCCTATACAGCTTGCTAATGATACAGTAAGTGTTCGGCCGGCGTTCACAACAAAGACCGCTATC ACCTAACACCACACCCGACACCCACTTAGGAGGGGAAAATTCGTCATCCATCAAGGCGATTAGGGTAGATTTACCAGGTGGGTTCAAGAC TTTAAGTTGCCATTTAGCTTGACCATCACACCGATCACCATTCGTTCTAAGCCTCTTAGTTAACCTGAGGCCGCTGTTTAGAAGACTGTA TCTTCAGTTGTAAACGTGGTGCTGGAGAGTTCCAAGTAGGTAATGGTGACTATCCGATCTTCAAATTGGGAGGGTAGTTAAGAGCACGAT CTCCGTCTCGTCGGTCCTGACTACGCGTTGCGTTTCCGCTCAGGCATTTTAGCCCGAGAGTGCGAGAAATTGCGCTACGCCCTTTATACT CTTGTGATTGTGGTAAACGCCACCGCAGTGAGGCCCCATGATGTCGCGCGGGTCTGTTTTCGTGGAGCTGATTTAAAGGGCCAACGATAC CAAGGGCAGTCTTTGTCATCATCGAGCGTGCGTAGGTAGGTGGCTGAAATTGGTAATAGCTGGCAGAACCTGACGATCCTAGAGAGCCGG GCTATTCCGACAGCCCAACACACACGCATACACTAGTGCGTTACGTCATTGGCGAAGGTGTACATTAGGCTTACTCCGTACGAAAAGGAT AGTCGGACTAAAGACATGTCGGCCGACCTCATGAATTCTATCCCCTGGTGCCGCCGTTCCCACGCTACAAGGTAATTACTTGGAGCGACA GCGACTCTGATCGCTCTTGTGCCTCAGAAGCAATCGAGAGACGATCCTCCTAGAACGAATATTGGCGGGCGACACAAGGATCGGGCTCAC AAGCTTGGTATCTGCCTCCGGGCTGTAGATGCAGCAGTACGTACAGTATCTGGCAAGTGGATAATGAGAATCAGCCTCGATATTACGTTT".split(' ')
k=len(pattern)
dis=0
for text in texts:
  hd=100000000
  i=0
  while i+k<=len(text):
    tp=text[i:i+k]
    c=0
    for j in range(k):
      if tp[j]!=pattern[j]:
        c=c+1
    if c<hd:
      hd=c
    i=i+1
  dis=dis+hd
print(dis)

# 3E: Contruct the De Bruijn Graph of a Collection of K Mers

def construct_de_bruijn_graph(kmers):
  graph={}
  for kmer in kmers:
    prefix=kmer[:-1]
    suffix=kmer[1:]
    if prefix not in graph:
      graph[prefix]=[]
    graph[prefix].append(suffix)
  return graph
#with open("input.txt","r") as file:
#  kmers=[line.strip() for line in file]
kmers=["GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"]
de_bruijn_graph=construct_de_bruijn_graph(kmers)
for node,neighbors in de_bruijn_graph.items():
  if neighbors:
    print(f"{node}->{','.join(neighbors)}")
  else:
    print(f"{node}->")

# 3G: Find an Eulerian Path in a Graph

from collections import defaultdict
def find_eulerian_path(graph):
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    for node in graph:
        out_degree[node] = len(graph[node])
        for neighbor in graph[node]:
            in_degree[neighbor] += 1
    start_node = None
    end_node = None
    for node in graph:
        if in_degree[node] + 1 == out_degree[node]:
            start_node = node
        if out_degree[node] + 1 == in_degree[node]:
            end_node = node
    if not start_node and not end_node:
        start_node = next(iter(graph))
    if end_node:
        graph[end_node].append(start_node)
        out_degree[start_node] += 1
    def visit(node):
        while graph[node]:
            visit(graph[node].pop())
        path.append(node)
    path = []
    visit(start_node)
    return path[::-1]
def read_graph_from_file(file_path):
    graph = defaultdict(list)
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(' -> ')
            node = int(parts[0])
            neighbors = parts[1].split(',')
            graph[node] = list(map(int, neighbors))
    return graph
file_path = "/content/rosalind_ba3g.txt"
graph = read_graph_from_file(file_path)
eulerian_path = find_eulerian_path(graph)
path_strings = []
for i in range(len(eulerian_path)):
    node = eulerian_path[i]
    if i == len(eulerian_path) - 1:
        path_strings.append(str(node))
    else:
        neighbors = ",".join(map(str, graph[node]))
        path_strings.append(f"{node}->{neighbors}")
output_line = "".join(path_strings)
print(output_line)

# 4B: Find Substrings of a Genome Encoding a Given Amino Acid String

#import Bio.Seq as seq
dna="ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
amino="MA"
k=len(amino)*3
for i in range(len(dna)-k+1):
  ts=seq.Seq(dna[i:i+k])
  tr=ts.reverse_complement()
  ta=ts.translate()
  tar=tr.translate()
  if ta==amino or tar==amino:
    print(ts)

# 4C: Generate Theoretical Spectrum of a Cyclic Peptide

pep="EAST"
spectrums=[pep,""]
amino_acid_masses={'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103,
                   'I': 113, 'L': 113, 'N': 114,'D': 115, 'K': 128, 'Q': 128, 'E': 129,
                   'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
n=len(pep)
pep=pep*2
for i in range(n):
  for j in range(n-1):
    spectrums.append(pep[i:i+j+1])
def getWeight(peptide):
  weight=0
  for acid in peptide:
    weight+=amino_acid_masses[acid]
  return weight
weights=[]
for spectrum in spectrums:
  weights.append(getWeight(spectrum))
weights.sort()
#for printing all the weights
# for weight in weights:
#   print(weight,end=" ")
#for printing the first 10 masses
for i in range(min(10,len(weights))):
  print(weights[i],end=" ")

# 4D: Compute the Number of Peptides of Given Total Mass

#import matplotlib.pyplot as plt
amino_acid_masses={'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103,
                   'I': 113, 'N': 114, 'D': 115, 'K': 128, 'E': 129, 'M': 131, 'H': 137,
                   'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def number_of_peptides(mass):
  dp=[0]*(mass+1)
  dp[0]=1
  for i in range(mass+1):
    for amin in amino_acid_masses:
      if i-amino_acid_masses[amin]>=0:
        dp[i]+=dp[i-amino_acid_masses[amin]]
  #return dp[mass]
  plt.figure(figsize=(20,10))
  plt.yscale('log')
  plt.plot(dp)
  plt.show()
number_of_peptides(1024)

# 4E: Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

import itertools
def cyclopeptide_sequencing(spec):
  mass=0
  result=[]
  for i in range(1,len(spec)):
    mass+=spec[i]
    result.append(spec[i])
    if mass==spec[len(spec)-1]:
      break
  result=list(itertools.permutations(result))
  return result
spec='''0 113 128 186 241 299 314 427'''.split()
spec=[eval(i) for i in spec]
res=cyclopeptide_sequencing(spec)
print(" ".join(map(str,res)))

# 4F: Compute the Score of a Cyclic Peptide Against a Spectrum

amino_acid_masses={'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103,
                   'I': 113, 'N': 114, 'D': 115, 'K': 128, 'E': 129, 'M': 131, 'H': 137,
                   'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def mass(peptide):
  return sum([amino_acid_masses[acid] for acid in peptide])
def cyclospectrum(peptide):
  double_peptide=peptide*2
  cspectrum=[0]
  for i in range(len(peptide)):
    for j in range(len(peptide)-1):
      cspectrum.append(mass(double_peptide[i:i+j+1]))
  cspectrum.append(mass(peptide))
  return cspectrum
def score(peptide,spectrum):
  actual_spectrum=cyclospectrum(peptide)
  s=0
  for i in range(len(spectrum)):
    for j in range(len(actual_spectrum)):
      if spectrum[i]==actual_spectrum[j]:
        s=s+1
        actual_spectrum[j]=-1
        break
  return s
peptide="EAST"
spectrum=[int(i) for i in "0 71 87 101 129 158 188 200 230 259 287 301 317 388".split()]
score(peptide,spectrum)

# 4H: Generate the Convolution of a Spectrum

def getconvol(cyclospectrum):
  difs={}
  for i in range(len(cyclospectrum)):
    for j in range(len(cyclospectrum)):
      d=cyclospectrum[i]-cyclospectrum[j]
      if d>0:
        if d in difs:
          difs[d]+=1
        else:
          difs[d]=1
  lis=[(difs[key],key) for key in difs]
  lis.sort()
  res=[]
  count_dict={}
  i=len(lis)-1
  while i>=0:
    k=lis[i][0]
    n=lis[i][1]
    for j in range(k):
      res.append(n)
      count_dict[n]=count_dict.get(n,0)+1
    i=i-1
  return res,count_dict
cyclospectrum=[int(i) for i in "0 71 87 101 129 158 188 200 230 259".split(' ')]
cnv,count_dict=getconvol(cyclospectrum)
print(" ".join(map(str,cnv)))
for key, count in count_dict.items():
    print(f"{key}: {count}")

# 5A: Find the Minimumm Number of Coins Needed to Make Change

def coinchange(coins,make):
  res=[0]*(make+1)
  for i in range(1,make+1):
    res[i]=min(res[i-coin]+1 for coin in coins if i>=coin)
  return res[make]
make=40
coins=1,5,10,20,25,50
coinchange(coins,make)

# 5C: Find a Longest Common Subsequence of Two Strings

s1,s2="AACCTTGG","ACACTGTGA"
n=len(s1)
m=len(s2)
dp=[[0]*(m+2) for i in range(n+1)]
for i in range(1,n+1):
  for j in range(1,m+1):
    if s1[i-1]==s2[j-1]:
      dp[i][j]=1+dp[i-1][j-1]
    else:
      dp[i][j]=max(dp[i-1][j],dp[i][j-1])
s=""
i=n
j=m
while i!=0 and j!=0:
  #print(i,j)
  if s1[i-1]==s2[j-1]:
    i=i-1
    j=j-1
    s=s+s1[i]
  elif dp[i][j]==dp[i-1][j]:
    i=i-1
  else:
    j=j-1
print(s[::-1])

# 5D: Find the Longest Path in a DAG

f=open("/content/rosalind_ba5d (2).txt")
src=int(f.readline())
sink=int(f.readline())
inf=-9999999
g=[[-9999999]*1001 for i in range(1001)]
n=0
for line in f:
  u,v=[int(i) for i in line.split(':')[0].split('->')]
  n=max(max(n,u),v)
  w=int(line.split(':')[1])
  g[u][v]=w
longest_path=[]
longest_length=inf
def dfs(cur,path,length):
  global longest_length
  global longest_path
  #print(path)
  if cur==sink:
    if length>longest_length:
      longest_length=length
      longest_path=path
    return
  for i in range(n+1):
    if g[cur][i]!=inf:
      dfs(i,path+[i],length+g[cur][i])
dfs(src,[0],0)
print(longest_length)
for i in longest_path:
  print(i,end="->")

# 5G: Compute the Edit Distance Between Two Strings

def edit_distance(X, Y):
  n = len(X)
  m = len(Y)
  dp = []
  for i in range(m+1):
    small = []
    for j in range(n+1):
      small.append(100000000)
    dp.append(small)
  dp[0][0] = 0
  n += 1
  m += 1
  for i in range(1, n):
    dp[0][i] = i
  for i in range(1, m):
    dp[i][0] = i
  X = "U" + X
  Y = "U" + Y
  for i in range(1, m):
    for j in range(1, n):
      if Y[i] == X[j]:
        dp[i][j] = dp[i-1][j-1]
      else:
        mn = min(dp[i][j-1], dp[i-1][j])
        dp[i][j] = min(dp[i-1][j-1], mn) + 1
  return dp[m-1][n-1]
if __name__ == "__main__":
  X = "PLEASANTLY"
  Y = "MEANLY"
  # with open('rosalind_ba5g.txt') as file:
  #   f = file.read().split()
  #   X = f[0]
  #   Y = f[1]
  d = edit_distance(X, Y)
  print(d)

# 6A: Implement Greedy Sorting to Sort a Permutation by Reversals

s=[int(i) for i in '-3 +4 +1 +5 -2'.split(' ')]
sorted=[abs(i) for i in s]
sorted.sort()
i=0
while i<len(s):
  j=i
  while(j<len(s) and abs(s[j])!=sorted[i]):
    j=j+1
  if i==j and s[i]==sorted[i]:
    i=i+1
    continue
  ts=[(-i) for i in s[i:j+1]]
  ts.reverse()
  s[i:j+1]=ts
  print('(',end="")
  for j in range(len(s)):
    if s[j]>=0:
      print('+',s[j],sep="",end="")
    else:
      print(s[j],end="")
    if j!=len(s)-1:
      print(end=" ")
  print(')')
  if s[i]==sorted[i]:
    i=i+1

# 6B: Compute the Number of Breakpoints in a Permutations

def number_of_breakpoints(p):
  p=[0]+p
  p.append(max(p)+1)
  num_bp=0
  for i in range(1,len(p)-1):
    if p[i]!=p[i-1]+1:
      num_bp+=1
  return num_bp
if __name__=="__main__":
  with open('/content/rosalind_ba6b.txt','r') as file:
    p=file.readline().strip()
    p=p.replace("(","").replace(")","")
    p=[int(x) for x in p.split()]
  print(number_of_breakpoints(p))

# 6E: Find All Shared K-Mers of a Pair of Strings

#import matplotlib.pyplot as plt
def reverse_comp(Seq):
  return Seq[::-1].translate(Seq.maketrans('ATCG','TAGC'))
def shared_kmers(k,seq1,seq2):
  result=[]
  seq1dict={seq1[i:i+k]: i for i in range(len(seq1)-k+1)}
  seq1dict.update({reverse_comp(seq1[i:i+k]): i for i in range(len(seq1)-k+1)})
  for j in range(len(seq2)-k+1):
    sub2=seq2[j:j+k]
    if sub2 in seq1dict:
      result.append([seq1dict[sub2],j])
    elif reverse_comp(sub2) in seq1dict:
      result.append([seq1dict[reverse_comp(sub2)],j])
  return result
if __name__=="__main__":
  with open('/content/rosalind_ba6e.txt','r') as file:
    input_lines=file.readlines()
  k=int(input_lines[0].rstrip())
  seq1=input_lines[1].rstrip()
  seq2=input_lines[2].rstrip()
  result=shared_kmers(k,seq1,seq2)
  for r in result:
    print('('+','.join(map(str,r))+')')

x_coordinates=[r[0] for r in result]
y_coordinates=[r[1] for r in result]
plt.scatter(x_coordinates,y_coordinates,marker='o',color='r')
plt.xlabel('position in sequence 1')
plt.ylabel('position in sequence 2')
plt.show()

# 6F: Implement Chromosome To Cycle

def chromo_to_cycle(chromosome):
  nodes=[]
  for block in chromosome:
    if block>0:
      nodes.append(2*block-1)
      nodes.append(2*block)
    else:
      nodes.append(-2*block)
      nodes.append(-2*block-1)
  return nodes
if __name__=="__main__":
  with open('/content/rosalind_ba6f.txt','r') as file:
    chromosome=file.readline().strip()
    chromosome=chromosome.replace("(","").replace(")","")
    chromosome=[int(x) for x in chromosome.split()]
  cycle=chromo_to_cycle(chromosome)
  print(" ".join(map(str,cycle)))

# 6G: Implement Cycle To Chromosome

def cycle_to_chromo(nodes):
  chromosome=[0]*(int(len(nodes)/2))
  for j in range(1,int(len(nodes)/2)+1):
    if nodes[2*j-1-1]<nodes[2*j-1]:
      chromosome[j-1]=nodes[2*j-1]/2
    else:
      chromosome[j-1]=-nodes[2*j-1-1]/2
  return [int(i) for i in chromosome]
s=[int(i) for i in '1 2 4 3 6 5 7 8'.split(' ')]
s=cycle_to_chromo(s)
print('(',end="")
for j in range(len(s)):
  if s[j]>=0:
    print('+',s[j],sep="",end="")
  else:
    print(s[j],end="")
  if j!=len(s)-1:
    print(end=" ")
print(')')

# 6H: Implement Colored Edges

import sys
def chromosome_to_cycle(Chromosome):
    Nodes = []
    for block in Chromosome:
        if block > 0:
            Nodes.append(2 * block - 1)
            Nodes.append(2 * block)
        else:
            Nodes.append(-2 * block)
            Nodes.append(-2 * block - 1)
    return Nodes
def colored_edges(P):
    Edges = list()
    for chromosome in P:
        Nodes = chromosome_to_cycle(chromosome)
        for j in range(1, len(Nodes), 2):
            if j != len(Nodes) - 1:
                Edges.append([Nodes[j], Nodes[j + 1]])
            else:
                Edges.append([Nodes[j], Nodes[0]])
    return Edges
if __name__ == "__main__":
    with open('/content/rosalind_ba6h.txt', 'r') as file:
        P = file.readline().strip()
        P = P[1:-1]
        P = P.split(')(')
        for i in range(len(P)):
            P[i] = [int(x) for x in P[i].split(' ')]
    result = colored_edges(P)
    for j in range(len(result)):
        result[j] = '(' + ', '.join(str(i) for i in result[j]) + ')'
    print(', '.join(result))

# 6J: Implement 2-BreakOnGenomeGraph

def TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4):
    if [i1, i2] in GenomeGraph:
        GenomeGraph.remove([i1, i2])
    else:
        GenomeGraph.remove([i2, i1])
    if [i3, i4] in GenomeGraph:
        GenomeGraph.remove([i3, i4])
    else:
        GenomeGraph.remove([i4, i3])
    GenomeGraph += [[i1, i3]] + [[i2, i4]]
    return GenomeGraph
if __name__ == "__main__":
    with open('/content/rosalind_ba6j.txt', 'r') as file:
        GenomeGraph = file.readline().strip()
        GenomeGraph = GenomeGraph[1:-1]
        GenomeGraph = GenomeGraph.split('), (')
        for i in range(len(GenomeGraph)):
            GenomeGraph[i] = GenomeGraph[i].split(', ')
            for j in range(len(GenomeGraph[i])):
                GenomeGraph[i][j] = int(GenomeGraph[i][j])
        i1, i2, i3, i4 = map(int, file.readline().strip().split(', '))
    result = TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4)
    for j in range(len(result)):
        result[j] = '(' + ', '.join(str(i) for i in result[j]) + ')'
    print(', '.join(result))