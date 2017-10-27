import random

mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}
rev = {}
for x in mapDict:
	if not mapDict[x] in rev:
		rev[mapDict[x]] = []
	rev[mapDict[x]].append(x)

class solution(object):
	def __init__(self,dna):
		self._dna=dna
		self._score=0
	def dna(self):
		return self._dna
	def score(self):
		return self._score
	def mutate(self,prob):
		for x in range(0,len(self._dna)):
			k = random.random()
			if k <= prob:
				self._dna[x] = random.choice(rev[self._dna[x]])
		return self
	def __eq__(self,other):
		return self.__dna == other._dna
				

def breed(a,b,ran_limit=100,prob=0.001):
	if not len(a.dna()) == len(b.dna()):
		return []
	children = []
	#half-half
	a1 = a.dna()[0:len(a.dna())/2]
	a2 = a.dna()[len(a.dna())/2:]
	b1 = b.dna()[0:len(b.dna())/2]
	b2 = b.dna()[len(b.dna())/2:]
	children = children + [solution(a1+b2),solution(b1+a2)]
	#zip
	a1 = []
	b1 = []
	for x in range(0,len(a.dna())):
		if x%2 == 0:
			a1.append(a.dna()[x])
			b1.append(b.dna()[x])
		else:
			a1.append(b.dna()[x])
			b1.append(b.dna()[x])
	children = children + [solution(a1),solution(b1)]
	#full random
	for x in range(0,ran_limit-4):
		a1 =[]
		for k in range(0,len(a.dna())):
			y = random.random()
			if y <=0.5:
				a1.append(a.dna()[k])
			else:
				a1.append(b.dna()[k])
		children.append(solution(a1))
	for x in children:
		x.mutate(prob)
	return children

def get_pool(a,size):
	children = []
	for x in range(0,size):
		a1 = []
		for k in range(0,len(a.dna())):
			a1.append(random.choice(rev[a.dna()[k]]))
		children.append(solution(a1))
	return children

def run(start,start_size=10,fitness,thresh,num_kids=10,mutate_prob=0.001):
	children = get_pool(start,start_size)
	for x in range(0,num_gen):
		next_gen = []
		for a in children:
			for b in children:
				if not a == b:
					next_gen = next_gen + breed(a,b,num_kids,mutate_prob)
		children = []
		for k in next_gen:
			if fitness(k) >= thresh:
				children.append(k)
	best = 0
	best_sol = None
	for k in children:
		fit = fitness(k)
		if fit > best:
			best = fit
			best_sol = k
	return (best_sol,best)
