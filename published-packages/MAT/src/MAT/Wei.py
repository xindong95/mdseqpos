#!/usr/bin/env python
# weili@jimmy.harvard.edu

import random, numpy, time, math
from __builtin__ import sum as suml


def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return numpy.array(result)

def sd(a):
    '''slow sd function
    '''
    alen = len(a) * 1.0
    amean = suml(a)/alen
    return math.sqrt(suml([(x-amean)*(x-amean) for x in a ])/(alen-1))

def variance(a):
    '''slow variance function
    '''
    alen = len(a) * 1.0
    amean = suml(a)/alen
    return suml([(x-amean)*(x-amean) for x in a ])/(alen-1)



def welch_ttest(a, b=''):
    '''welch_test(a, b): return Welch  t-test statistics
    '''
    alen = len(a)
    avar = variance(a)
    if len(b):
        blen = len(b)
        bvar = variance(b)
        t = (suml(a)/alen - suml(b)/blen)/math.sqrt(avar/alen + bvar/blen)
        #df = math.pow((avar/alen + bvar/blen), 2) / (math.pow(avar/alen, 2)/(alen-1) + math.pow(bvar/blen, 2)/(blen-1) )
    else:       # one sample
        t = suml(a)/alen/math.sqrt(avar/alen)
    return t




RCLookup = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'A', 'N':'N','-':'-',
            'a':'T', 'c':'G', 'g':'C', 't':'A', 'u':'A', 'n':'N'}
def RevComp(seqs):
    try:
        rseq = map(lambda x: (RCLookup[x]), list(seqs))
        rseq.reverse()

        return "".join(rseq)
    except:
        print "error"
        return "";


class Stack(object):
    '''stack object
    Stack(start = [])
    '''
    def __init__(self, start=[]):            # self is the instance object
        self.stack = []                      # start is any sequence: stack..
        for x in start: self.push(x)
    def push(self, obj):                     # methods: like module + self
        self.stack.append(obj)               # top is end of list
    def pop(self):
        if not self.stack: raise error, 'underflow'
        return self.stack.pop(  )              # like fetch and delete stack[-1]
    def top(self):
        if not self.stack: raise error, 'underflow'
        return self.stack[-1]
    def empty(self):
        return not self.stack                # instance.empty(  )
    def __getitem__(self, offset):
        return self.stack[offset]            # instance[offset], in, for

    def __repr__(self):
        return '[Stack:%s]' % self.stack          # print, backquotes,..
    def __cmp__(self, other):
        return cmp(self.stack, other.stack)       # '==', '>, '<=', '!=',..
    def __len__(self):
        return len(self.stack)                    # len(instance), not instance
    def __add__(self, other):
        return Stack(self.stack + other.stack)    # instance1 + instance2
    def __mul__(self, reps):
        return Stack(self.stack * reps)           # instance * reps
    def __getitem__(self, offset):
        return self.stack[offset]                 # intance[offset], in, for
    def __getslice__(self, low, high):
        return Stack(self.stack[low : high])      # instance[low:high]
    def __getattr__(self, name):
        return getattr(self.stack, name)          # instance.sort()/reverse(  )/..


class Set(object):
    def __init__(self, value = []):
        self.data = {}                     # manages a local dictionary
        self.concat(value)                 # hashing: linear search times
    def intersect(self, other):
        res = {}
        for x in other:                    # other: a sequence or Set
            if self.data.has_key(x):       # use hash-table lookup
                res[x] = None
        return Set(res.keys())             # a new dictionary-based Set
    def union(self, other):
        res = {}                           # other: a sequence or Set
        for x in other:                    # scan each set just once
            res[x] = None
        for x in self.data.keys():         # '&' and '|' come back here
            res[x] = None                  # so they make new fastset's
        return Set(res.keys())
    def concat(self, value):
        for x in value: self.data[x] = None

    # inherit and, or, len
    def __getitem__(self, key):  return self.data.keys()[key]
    def __repr__(self):          return '<Set:' + `self.data.keys()` + '>'
    def __len__(self):          return len(self.data)
    def __and__(self, other):   return self.intersect(other)
    def __or__(self, other):    return self.union(other)




class Graph(object):
    '''Node in graph.
    A = Graph('A')
    B= Graph('B')
    C = Graph('C')
    A.arcs = [B, C]
    A.search(C)
    '''
    def __init__(self, label, extra=None):
        self.name = label                                # nodes=inst objects
        self.data = extra                                # graph=linked objs
        self.arcs = []
    def __repr__(self):
        return self.name
    def search(self, goal):
        Graph.solns = []
        self.generate([self], goal)
        Graph.solns.sort(lambda x,y: cmp(len(x), len(y)))
        return Graph.solns
    def generate(self, path, goal):
        if self == goal:                                 # class == tests addr
            Graph.solns.append(path)                     # or self.solns: same
        else:
            for arc in self.arcs:
                if arc not in path:
                    arc.generate(path + [arc], goal)


'''
xpermutations takes all elements from the sequence, order matters.

xcombinations takes n distinct elements from the sequence, order matters.

xuniqueCombinations takes n distinct elements from the sequence, order is irrelevant.

xselections takes n elements (not necessarily distinct) from the sequence, order matters.

'''

def xcombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcombinations(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc

def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc

def xselections(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss

def xpermutations(items):
    return xcombinations(items, len(items))






