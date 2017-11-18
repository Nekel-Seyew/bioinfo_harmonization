import sys
import random

codons = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}
aa = {}
for c in codons:
    if codons[c] in aa:
        aa[codons[c]].append(c)
    else:
        aa[codons[c]] = [c]
host_freq = sys.argv[1]
orig_freq = sys.argv[2]
seq_file = sys.argv[3]

host_dict = {}
orig_dict = {}
host_rel_freq = {}
orig_rel_freq = {}
with open(host_freq) as f:
	for line in f:
		line = line.rstrip().rsplit()
		host_dict[line[0]] = float(line[1])
        #host_dict[line[0]] = random.randrange(1, 100)

for a in aa:
    ab_freq = 0
    for c in aa[a]:
        ab_freq += (host_dict[c] / 10)
    factor = 100 / ab_freq
    for c in aa[a]:
        host_rel_freq[c] = host_dict[c] * factor / 10
    print ab_freq
    print factor
print host_rel_freq

with open(orig_freq) as f:
	for line in f:
		line = line.rstrip().rsplit()
		orig_dict[line[0]] = float(line[1])
        #orig_dict[line[0]] = random.randrange(1, 100)

for a in aa:
    ab_freq = 0
    for c in aa[a]:
        ab_freq += (orig_dict[c] / 10)
    factor = 100 / ab_freq
    for c in aa[a]:
        orig_rel_freq[c] = orig_dict[c] * factor / 10
    print ab_freq
    print factor
print orig_rel_freq

seq = ""
with open(seq_file) as f:
    for line in f:
        if ">" not in line:
            line = line.rstrip()
            seq += line
print seq
print len(seq)
def find_cutoff():
    n_clust = int(len(seq) / 90)
    print n_clust

    clust_array = []
    clusts = 0
    cutoff = 10
    cutoffs = []
    while cutoff < 21:
        clust_array = []
        for i in range(0, int(len(seq) / 3) * 3, 3):
            #print i/3, seq[i:i+3], orig_rel_freq[seq[i:i+3]]
            if orig_rel_freq[seq[i:i+3]] <= cutoff:
                clust_array.append(i / 3)
        n = 1
        for i in range(1, len(clust_array) - 1):
            if i == 1 and clust_array[i - 1] < clust_array[i] - 14:
                n -= 1
            if clust_array[i] > clust_array[i-1] + 14 and clust_array[i] > clust_array[i + 1] - 14:
                n += 1
        cutoffs.append((cutoff, n))
        cutoff += 0.5
    print cutoffs
    cutoffs.sort(key=lambda x:abs(x[1] - n_clust))
    return cutoffs[0][0]

x = find_cutoff()
print x
