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
	def give_score(self,score):
		self._score = score
	def score(self):
		return self._score
	def mutate(self,prob):
		for x in range(0,len(self._dna)):
			k = random.random()
			if k <= prob:
				self._dna[x] = random.choice(rev[mapDict[self._dna[x]]])
		return self
	def mutate_pos(self,pos):
		self._dna[pos] = random.choice(rev[mapDict[self._dna[pos]]])
		return self
	def __eq__(self,other):
		return self._dna == other._dna
	def __repr__(self):
		return " Score: "+str(self._score)
				

def breed(a,b,ran_limit=100,prob=0.001):
	if not len(a.dna()) == len(b.dna()):
		return []
	children = []
	#half-half
	a1 = a.dna()[0:int(len(a.dna())/2)]
	a2 = a.dna()[int(len(a.dna())/2):]
	b1 = b.dna()[0:int(len(b.dna())/2)]
	b2 = b.dna()[int(len(b.dna())/2):]
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
	#two at a time
	a1=[]
	b1=[]
	switch=False
	for x in range(0,len(a.dna())):
		if x%2 == 0:
			switch = not switch
		if switch:
			a1.append(a.dna()[x])
			b1.append(b.dna()[x])
		else:
			a1.append(b.dna()[x])
			b1.append(a.dna()[x])
	children = children + [solution(a1),solution(b1)]
	#three at a time
	a1=[]
	b1=[]
	switch=False
	for x in range(0,len(a.dna())):
		if x%3 == 0:
			switch = not switch
		if switch:
			a1.append(a.dna()[x])
			b1.append(b.dna()[x])
		else:
			a1.append(b.dna()[x])
			b1.append(a.dna()[x])
	children = children + [solution(a1),solution(b1)]
	#10 at a time
	a1=[]
	b1=[]
	switch=False
	for x in range(0,len(a.dna())):
		if x%10 == 0:
			switch = not switch
		if switch:
			a1.append(a.dna()[x])
			b1.append(b.dna()[x])
		else:
			a1.append(b.dna()[x])
			b1.append(a.dna()[x])
	children = children + [solution(a1),solution(b1)]

	#full random
	for x in range(0,ran_limit-10):
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
			a1.append(random.choice(rev[mapDict[a.dna()[k]]]))
		children.append(solution(a1))
	return children

def avg(a):
	k = 0;
	for i in a:
		k +=i
	return k/len(a)

#fitness function must return a number, and lower is better
def run(start,fitness,worst_genes,start_size=10,num_kids=10,num_gen=100,mutate_prob=0.025):
	print("Creating first gen pool")
	children = get_pool(start,start_size)
	print("Starting generations")
	last_best = fitness(start)
	window = []
	best = 2**32
	best_sol = None
	for x in range(0,num_gen):
		print("On generation: ",x)
		next_gen = []
		for a in children:
			for b in children:
				if not a == b:
					next_gen = next_gen + breed(a,b,num_kids,mutate_prob)
		if len(next_gen) == 0:
			for k in children:
				k.mutate(mutate_prob)
			continue
		for a in next_gen:
			a.give_score(fitness(a))
		next_gen = sorted(next_gen,key=lambda x: x.score())
		children = []
		i = 0
		for a in next_gen:
			if not a in children and i < start_size:
				children.append(a)
				i += 1
		if children[0].score() < best:
			best = children[0].score()
			best_sol = children[0]
		print("the children: ",children)
		if len(window) < 10:
			window.append(children[0].score())
		else:
			window = window[1:] + [children[0].score()]
		if abs(avg(window) - children[0].score()) < 2:
			for k in children:
				worst = worst_genes(k)
				for i in worst:
					k.mutate_pos(i)
	return (best_sol.dna(),best)
