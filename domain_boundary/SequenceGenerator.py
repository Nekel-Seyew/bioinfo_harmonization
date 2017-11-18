import sys
class SequenceGenerator:
    def __init__(self, host, orig, seq):
        self.host_freq = host
        self.orig_freq = orig
        self.seq_file = seq

        self.host_dict = {}
        self.orig_dict = {}
        self.host_rel_freq = {}
        self.orig_rel_freq = {}

        self.seq = ""
        self.new_seq = ""
        self.rare_cutoff = -1

        self.clustered_site = []
        self.isolated_site = []

        self.codons = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}

        self.end_link_aas = ["Y", "H", "W", "I", "L", "V", "S", "T", "P", "C"]

        self.aa = {}
        for c in self.codons:
            if self.codons[c] in self.aa:
                self.aa[self.codons[c]].append(c)
            else:
                self.aa[self.codons[c]] = [c]

    def load_files(self):
        with open(self.host_freq) as f:
            for line in f:
                line = line.rstrip().rsplit()
                self.host_dict[line[0]] = float(line[1])
                #host_dict[line[0]] = random.randrange(1, 100)

        for a in self.aa:
            ab_freq = 0
            for c in self.aa[a]:
                ab_freq += (self.host_dict[c] / 10)
            factor = 100 / ab_freq
            for c in self.aa[a]:
                self.host_rel_freq[c] = self.host_dict[c] * factor / 10

        with open(self.orig_freq) as f:
            for line in f:
                line = line.rstrip().rsplit()
                self.orig_dict[line[0]] = float(line[1])
                #orig_dict[line[0]] = random.randrange(1, 100)

        for a in self.aa:
            ab_freq = 0
            for c in self.aa[a]:
                ab_freq += (self.orig_dict[c] / 10)
            factor = 100 / ab_freq
            for c in self.aa[a]:
                self.orig_rel_freq[c] = self.orig_dict[c] * factor / 10

        with open(self.seq_file) as f:
            for line in f:
                if ">" not in line:
                    line = line.rstrip()
                    self.seq += line
        print "Original Sequence in " + self.orig_freq[:-9]
        print self.seq

    def find_cutoff(self):
        n_clust = int(len(self.seq) / 90)

        clust_array = []
        clusts = 0
        cutoff = 10
        cutoffs = []
        while cutoff < 21:
            clust_array = []
            for i in range(0, int(len(self.seq) / 3) * 3, 3):
                #print i/3, seq[i:i+3], orig_rel_freq[seq[i:i+3]]
                if self.orig_rel_freq[self.seq[i:i+3]] <= cutoff:
                    clust_array.append(i / 3)
            n = 1
            for i in range(1, len(clust_array) - 1):
                if i == 1 and clust_array[i - 1] < clust_array[i] - 14:
                    n -= 1
                if clust_array[i] > clust_array[i-1] + 14 and clust_array[i] > clust_array[i + 1] - 14:
                    n += 1
            cutoffs.append((cutoff, n))
            cutoff += 0.5
        #print cutoffs
        cutoffs.sort(key=lambda x:abs(x[1] - n_clust))
        self.rare_cutoff = cutoffs[0][0]

        clust_array = []
        for i in range(0, int(len(self.seq) / 3) * 3, 3):
            #print i/3, seq[i:i+3], orig_rel_freq[seq[i:i+3]]
            if self.orig_rel_freq[self.seq[i:i+3]] <= self.rare_cutoff:
                clust_array.append(i / 3)

        for i in range(0, len(clust_array)):
            if i == 0:
                if clust_array[i+1] - 14 > clust_array[i]:
                    self.isolated_site.append(clust_array[i] * 3)
                else:
                    self.clustered_site.append(clust_array[i] * 3)
            elif i > 0 and i < len(clust_array) - 1:
                if clust_array[i-1] + 14 < clust_array[i] and clust_array[i+1] - 14 > clust_array[i]:
                    self.isolated_site.append(clust_array[i] * 3)
                else:
                    self.clustered_site.append(clust_array[i] * 3)
            else:
                if clust_array[i-1] + 14 < clust_array[i]:
                    self.isolated_site.append(clust_array[i] * 3)
                else:
                    self.clustered_site.append(clust_array[i] * 3)

    def find_closest_lower(self, c, freq):
        a = self.codons[c]
        candidates = self.aa[a]
        cand_freqs = []
        for cand in candidates:
            cand_freqs.append((cand, self.host_rel_freq[cand]))
        cand_freqs.sort(key = lambda x:x[1])
        #print "Codon: " + c + "\tFreq: " + str(freq)
        #print cand_freqs

        i = 0
        new_freq = cand_freqs[i][1]
        while new_freq <= freq:
            i += 1
            new_freq = cand_freqs[i][1]
        
        #print cand_freqs[i-1][0]
        return cand_freqs[i-1][0]

    def find_closest(self, c, freq):
        a = self.codons[c]
        candidates = self.aa[a]
        cand_freqs = []
        for cand in candidates:
            cand_freqs.append((cand, self.host_rel_freq[cand]))
        cand_freqs.sort(key = lambda x:x[1])
        #print "Codon: " + c + "\tFreq: " + str(freq)
        #print cand_freqs
        i = 0
        new_freq = cand_freqs[i][1]
        while new_freq < freq and i < len(cand_freqs) - 1:
            i += 1
            new_freq = cand_freqs[i][1]

        if abs(freq - cand_freqs[i-1][1]) < abs(freq - cand_freqs[i][1]):
            #print cand_freqs[i-1][0]
            return cand_freqs[i-1][0]
        else:
            #print cand_freqs[i][0]
            return cand_freqs[i][0]

    def replace_codons(self):
        for i in range(0, int(len(self.seq)/3) * 3, 3):
            c = self.seq[i:i+3]
            aa = self.codons[c]
            freq = self.orig_rel_freq[c]
            if freq <= self.rare_cutoff and aa in self.end_link_aas and i in self.clustered_site:
                new_c = self.find_closest_lower(c, freq)
                self.new_seq += new_c
            else:
                new_c = self.find_closest(c, freq)
                self.new_seq += new_c

    def run(self):
        self.load_files()
        self.find_cutoff()
        self.replace_codons()
        return self.new_seq

if __name__ == '__main__':
    x = SequenceGenerator(sys.argv[1], sys.argv[2], sys.argv[3])
    x.run()
    print "\nSequence harmonized to " + x.host_freq[:-9]
    print x.new_seq
