import random
import networkx as nx
import time
import math

from multiprocessing import Pool
import dill as pickle

import sys

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

primes =[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199]

def factors(num):
    i = 0
    a = []
    while num > 1:
        if num % primes[i] == 0:
            num /= primes[i]
            a.append(primes[i])
        else:
            i += 1
    return a

def mulsum(lst):
    i = 1
    for x in lst:
        i *= x
    return i

class solution:
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
    	return tuple(self._dna) == tuple(other._dna)
    def __repr__(self):
    	return " Score: "+str(self._score)
    def __hash__(self):
                return hash(tuple(self._dna))
    			
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

def getNeighbors(graph,ni):
    a = []
    for k in ni:
        a.append(graph.node[k]['sol'])
    return a
def getSols(graph):
    a = []
    for x in graph.nodes():
        a.append(graph.node[x]['sol'])
    return a

def get_best_neighbor(graph,ni):
    best = 10000000000000000
    best_sol = None
    for k in ni:
        if graph.node[k]['sol'].score() < best:
            best = graph.node[k]['sol'].score()
            best_sol = graph.node[k]['sol']
    return best_sol

def get_rand_neighbor(graph,ni):
    neigh = []
    for k in ni:
        neigh.append(graph.node[k])
    ret = None if len(neigh) == 0 else random.choice(neigh)['sol']
    return ret

def multiproc_vertex(tupin):
    graph = pickle.loads(tupin[0])
    ve = tupin[1]
    fitness = pickle.loads(tupin[2])
    num_kids = tupin[3]
    mutate_prob = tupin[4]
    kids = []
    neighbors_index = graph.neighbors(ve)
    neighbors = getNeighbors(graph,neighbors_index)
    for ne in neighbors:
        #print(graph.node[ne]['sol'])
        kids = kids + breed(graph.node[ve]['sol'],ne,num_kids,mutate_prob)
    #rand_neighbor = get_rand_neighbor(graph,neighbors_index)
    #kids = [] if rand_neighbor is None else breed(rand_neighbor,graph.node[ve]['sol'],num_kids,mutate_prob)
    for ki in kids:
        ki.give_score(fitness(ki))
    kids = sorted(kids,key=lambda x: x.score())
    #print(kids)
    if len(kids) > 0 and kids[0].score() < graph.node[ve]['sol'].score():
        return (ve,kids[0])


def graph_run(start,fitness,worst_genes,start_size=10,num_kids=10,num_gen=100,mutate_prob=0.025,verts=20):
    starttime = time.perf_counter()
    children = get_pool(start,verts)
    for a in children:
        a.give_score(fitness(a))
    children = sorted(children, key=lambda x: x.score())
    #now we need to make the graph
    lst = factors(verts)
    #graph = nx.complete_graph(verts)
    #graph = nx.desargues_graph()
    #graph = nx.dodecahedral_graph()
    half = math.ceil(len(lst)/2)
    #graph = nx.watts_strogatz_graph(verts,10,0.1)
    #graph = nx.windmill_graph(mulsum(lst[half:]),mulsum(lst[:half]))
    graph = nx.windmill_graph(mulsum(lst[:half]),mulsum(lst[half:]))
    #graph = nx.caveman_graph(mulsum(lst[:half]),mulsum(lst[half:]))
    #graph = nx.caveman_graph(mulsum(lst[half:]),mulsum(lst[:half]))
    #graph = nx.gnp_random_graph(verts,0.1)
    #graph = nx.erdos_renyi_graph(verts,0.1)
    #graph = nx.barabasi_albert_graph(verts,int(verts/3))
    #graph = nx.hypercube_graph(6)
    #graph = nx.grid_2d_graph(mulsum(lst[:half]),mulsum(lst[half:]))
    graph = nx.convert_node_labels_to_integers(graph)
    gattr = {}
    #print(graph.nodes())
    for n in graph.nodes():
        gattr[n] = {'sol':children[n]}
        #graph.add_node(n,sol=children[n])
    nx.set_node_attributes(graph,gattr)
    #graph now done, sheeesh
    for x in range(num_gen):
        #replaces=[]
        nodetupes = []
        for ve in graph.nodes():
            nodetupes += [(pickle.dumps(graph),ve,pickle.dumps(fitness),num_kids,mutate_prob)]
            #kids = []
            #neighbors_index = graph.neighbors(ve)
            #neighbors = getNeighbors(graph,neighbors_index)
            #for ne in neighbors_index:
            #    kids = kids + breed(graph.node[ve]['sol'],graph.node[ne]['sol'],num_kids,mutate_prob)
            #rand_neighbor = get_rand_neighbor(graph,neighbors_index)
            #kids = [] if rand_neighbor is None else breed(rand_neighbor,graph.node[ve]['sol'],num_kids,mutate_prob)
            #for ki in kids:
            #    ki.give_score(fitness(ki))
            #kids = sorted(kids,key=lambda x: x.score())
            #if len(kids) > 0 and kids[0].score() < graph.node[ve]['sol'].score():
            #    replaces.append((ve,kids[0]))
        graphsols = set()
        with Pool(16) as p:
            replaces = p.map(multiproc_vertex,nodetupes)
        for ve in graph.nodes():
            graphsols.add(graph.node[ve]['sol'])
        for rep in replaces:
            if rep == None:
                continue
            insol = False
            insol = rep[1] in graphsols
            #for ve in graph.nodes():
            #    if graph.node[ve]['sol'] == rep[1]:
            #        insol = True
            if not insol:
                graph.node[rep[0]]['sol']=rep[1]
        #print("GEN %i "%x,getSols(graph))
        #nx.set_node_attributes(graph,gattr)
    best_sol = sorted(getSols(graph),key=lambda x: x.score())[0]
    endtime = time.perf_counter()
    #print(endtime-starttime)
    return (best_sol.dna(),best_sol.score(),endtime-starttime)



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
                #breeding stage
                
    	for a in children:
    		for b in children:
    			if not a == b:
    				next_gen = next_gen + breed(a,b,num_kids,mutate_prob)
    	if len(next_gen) == 0:
    		for k in children:
    			k.mutate(mutate_prob)
    		continue
                #scoring stage
    	for a in next_gen:
    		a.give_score(fitness(a))
    	next_gen = sorted(next_gen,key=lambda x: x.score())
    	children = []
    	i = 0
                #reduce number of children
    	for a in next_gen:
    		if not a in children and i < start_size:
    			children.append(a)
    			i += 1
    	if len(children) > 0 and children[0].score() < best:
    		best = children[0].score()
    		best_sol = children[0]
                #print children
    	print("the children: ",children)
    	if len(window) < 10:
    		window.append(children[0].score())
    	else:
    		window = window[1:] + [children[0].score()]
                #if not making progress, do something about it
    	if abs(avg(window) - children[0].score()) < 2:
    		for k in children:
    			worst = worst_genes(k)
    			for i in worst:
    				k.mutate_pos(i)
    return (best_sol.dna(),best)
