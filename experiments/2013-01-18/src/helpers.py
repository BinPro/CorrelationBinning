#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint 
from probin.model.composition import multinomial as mn

def sample_contigs(genome, n, min_length, max_length):
    contigs = []
    identifier = ">" + genome.id
    for i in range(n):
        l = randint(min_length, max_length)
        start = randint(0, (len(genome)-l))
        contig = Seq(genome.seq[start:start+l])
        contigs.append(identifier + " contig_number " + str(i) + "\n" + str(genome.seq[start:start+l]))
    return contigs

def sample_contig(genome, min_length, max_length):
    l = randint(min_length, max_length)
    start = randint(0, (len(genome)-l))
    contig = SeqRecord(genome.seq[start:start+l])
    rest = SeqRecord(genome.seq[0:start] + genome.seq[start+l:-1])
    return contig, rest

def score_contig(signature, gen_par, group_pars, rest_pars):
    gen_score = mn.log_probability(signature, gen_par)
    group_scores = []
    rest_scores = []
    for par in group_pars:
        group_scores.append(mn.log_probability(signature, par))
    for par in rest_pars:
        rest_scores.append(mn.log_probability(signature,par))
    return gen_score, group_scores, rest_scores
