#! /usr/bin/env python

import numpy

def sample(p):
    """Given a NumPy array p, select an index based on the probabilities in p."""
    assert type(p) is numpy.ndarray
    assert p.ndim == 1
    assert p.shape[0] > 0
    assert (p >= 0).sum() == p.shape[0]
    if p.sum() > 0:
        r = numpy.random.random()
        c = p.cumsum() / p.sum()
        i = numpy.where(r < c)[0].min()
    else:
        i = numpy.random.randint(p.shape[0])
    return i

def read_freq(motif, base=('A', 'C', 'G', 'T')):
    """Read frequencies from the output of MDscan."""
    freq = motif
    return freq, base

def fasta(motif, fastafilename, n=10000, base=('A', 'C', 'G', 'T'), verbose=True):
    """Convert output of MDscan to a list of sequences with the correct frequencies."""
    freq, base = read_freq(motif, base)
    freq = numpy.array(freq)
    fastafile = open(fastafilename, 'w')
    for seqnum in xrange(n):
        seqstr = "".join(base[sample(p)] for p in freq)
        print >> fastafile, ">%d" % seqnum
        print >> fastafile, seqstr
    fastafile.close()

