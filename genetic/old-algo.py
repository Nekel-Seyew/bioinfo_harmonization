import random
import networkx as nx

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
    	return tuple(self._dna) == tuple(other._dna)
    def __repr__(self):
    	return " Score: "+str(self._score)
    def __hash__(self):
                return hash(tuple(self._dna))
    			

class graph(object):
        def __init__(self,vert_num):
            self._graph = [[0 for y in range(vert_num)] for x in range(vert_num)]
            self._vert_map = {}
            self._verts = 0
            self._inv_vert_map = {}

        def num_vert(self):
            return len(self._graph)
        def give_vert_obj(self,obj):
            self._inv_vert_map[self._verts] = obj
            self._vert_map[obj] = self._verts
            self._verts += 1
        def give_edge(self,a,b):
            if type(a) is int and type(b) is int:
                self._graph[a][b] = 1
                self._graph[b][a] = 1
            else:
                if not a in self._vert_map:
                    print("ERROR, ",a," is not in the graph!")
                if not b in self._vert_map:
                    print("ERROR, ",b," is not in the graph!")
                ai = self._vert_map[a]
                bi = self._vert_map[b]
                self._graph[ai][bi] = 1
                self._graph[bi][ai] = 1
        def get_neighbors(self,a):
            if not a in self._vert_map:
                print("ERROR, ",a," IS NOT IN THE GRAPH!")
            ret = []
            ai = a if type(a) is int else self._vert_map[a]
            for k in range(len(self._graph)):
                if self._graph[ai][k] == 1:
                    ret.append(self._inv_vert_map[k])
            return ret
        def replace_vert(self,old,new):
            if new in self._vert_map:
                return False
            ai = self._vert_map[old]
            del self._vert_map[old]
            self._vert_map[new] = ai
            self._inv_vert_map[ai] = new
        def get_verts(self):
            return [key for key in self._vert_map]


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



def graph_run(start,fitness,worst_genes,start_size=10,num_kids=10,num_gen=100,mutate_prob=0.025,verts=20):
    children = get_pool(start,verts)
    for a in children:
        a.give_score(fitness(a))
    print("Sorting children, there are: ",len(children))
    children = sorted(children, key=lambda x: x.score())[:verts]
    #now we need to make the graph
    a = children[0]
    b = children[1]
    c = children[2]
    d = children[3]
    e = children[4]
    f = children[5]
    g = children[6]
    h = children[7]
    i = children[8]
    j = children[9]
    k = children[10]
    l = children[11]
    m = children[12]
    n = children[13]
    o = children[14]
    p = children[15]
    q = children[16]
    r = children[17]
    s = children[18]
    t = children[19]
    le_graph = graph(verts)
    le_graph.give_vert_obj(a)
    le_graph.give_vert_obj(b)
    le_graph.give_vert_obj(c)
    le_graph.give_vert_obj(d)
    le_graph.give_vert_obj(e)
    le_graph.give_vert_obj(f)
    le_graph.give_vert_obj(g)
    le_graph.give_vert_obj(h)
    le_graph.give_vert_obj(i)
    le_graph.give_vert_obj(j)
    le_graph.give_vert_obj(k)
    le_graph.give_vert_obj(l)
    le_graph.give_vert_obj(m)
    le_graph.give_vert_obj(n)
    le_graph.give_vert_obj(o)
    le_graph.give_vert_obj(p)
    le_graph.give_vert_obj(q)
    le_graph.give_vert_obj(r)
    le_graph.give_vert_obj(s)
    le_graph.give_vert_obj(t)
    le_graph.give_edge(a,b)
    le_graph.give_edge(a,j)
    le_graph.give_edge(a,k)
    le_graph.give_edge(b,c)
    le_graph.give_edge(b,l)
    le_graph.give_edge(c,d)
    le_graph.give_edge(c,m)
    le_graph.give_edge(d,e)
    le_graph.give_edge(d,n)
    le_graph.give_edge(e,o)
    le_graph.give_edge(e,f)
    le_graph.give_edge(f,p)
    le_graph.give_edge(f,g)
    le_graph.give_edge(g,q)
    le_graph.give_edge(g,h)
    le_graph.give_edge(h,r)
    le_graph.give_edge(h,i)
    le_graph.give_edge(i,s)
    le_graph.give_edge(i,j)
    le_graph.give_edge(j,t)
    le_graph.give_edge(k,l)
    le_graph.give_edge(l,m)
    le_graph.give_edge(m,n)
    le_graph.give_edge(n,o)
    le_graph.give_edge(o,p)
    le_graph.give_edge(p,q)
    le_graph.give_edge(q,r)
    le_graph.give_edge(r,s)
    le_graph.give_edge(s,t)
    le_graph.give_edge(t,k)
    #graph now done, sheeesh
    for x in range(num_gen):
        verts = le_graph.get_verts()
        replaces=[]
        for ve in verts:
            kids = []
            neighbors = le_graph.get_neighbors(ve)
            for ne in neighbors:
                kids = kids + breed(ve,ne,num_kids,mutate_prob)
            for ki in kids:
                ki.give_score(fitness(ki))
            kids = sorted(kids,key=lambda x: x.score())
            if kids[0].score() < ve.score():
                replaces.append((ve,kids[0]))
        for rep in replaces:
            le_graph.replace_vert(rep[0],rep[1])
        print("GEN ",x," Current verts: ",le_graph.get_verts())
    best_sol = sorted(le_graph.get_verts(),key=lambda x: x.score())[0]
    return (best_sol.dna(),best_sol.score())



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
    	if children[0].score() < best:
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
