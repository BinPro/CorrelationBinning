#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint 
from probin.model.composition import multinomial as mn
from probin.dna import DNA

def sample_contigs(genome, n, min_length, max_length):
    contigs = []
    identifier = ">" + genome.id
    for i in range(n):
        l = randint(min_length, max_length)
        start = randint(0, (len(genome)-l))
        contig = Seq(genome.seq[start:start+l])
        contigs.append(identifier + " contig_number " + str(i) + "\n" + str(genome.seq[start:start+l]))
    return contigs

def sample_contig(genome, x_st):
    min_length = x_st.contig_min_length
    max_length = x_st.contig_max_length
    l = randint(min_length, max_length)
    start = randint(0, (len(genome.seq)-l))
    contig = DNA(id = genome.id + " contig", seq = genome.seq[start:start+l])
    rest = DNA(id= genome.id, seq = genome.seq[0:start] + genome.seq[start+l:-1])
    return contig, rest

def score_contig(contig, gen, genome_index, test):
    score_obj = ScoreCollection()
    score_gen = Score(mn.log_probability(contig.signature, gen.par()), gen.id)
    score_obj.genome = score_gen
    group_genomes = test.group.all_genomes_but_index(genome_index)
    for gen in group_genomes:
        score_obj.group.append(Score(mn.log_probability(contig.signature, gen.par()),gen.id)) 
    for group in test.rest_groups:
        outside_group = []
        for gen in group.genomes:
            score = Score(mn.log_probability(contig.signature,gen.par()), gen.id)
            outside_group.append(score)
        score_obj.other.append(outside_group)
    return score_obj

def all_but_index(l,i):
    return l[0:i] + l[i+1:]

class GenomeGroup:
    """A wrapper for the genome groups used in pairwise testing"""
    def __init__(self, name):
        self.genomes = []
        self.name = name
    def add_genome(self, genome):
        self.genomes.append(genome)
    def parameters(self):
        return [genome.par for genome in self.genomes]
    def __len__(self):
        return len(self.genomes)
    def all_genomes_but_index(self,i):
        return self.genomes[0:i] + self.genomes[i+1:]


def par(self):
    try:
        return self.parameters
    except:
        pa = {}
        n = sum(self.signature.values())
        for i,cnt in self.signature.items():
            pa[i] = cnt/float(n)
        self.parameters = pa
        return pa

DNA.par = par

class ExperimentSetting:
    """A class for storing variables related to all experiments"""
    def __init__(self, prio, mode, no_contigs, contig_min_length, contig_max_length):
        self.prio = prio
        self.mode = mode
        self.no_contigs = no_contigs
        self.contig_min_length = contig_min_length
        self.contig_max_length = contig_max_length

class Test:
    def __init__(self,x_st, group, rest_groups):
        self.n = len(group)
        self.count_per_g = self.count_per_g(x_st.no_contigs, self.n)
        self.group = group
        self.rest_groups = rest_groups
        self.x_st = x_st

    def count_per_g(self, no_contigs, n):
        return [round(no_contigs/float(self.n))]*self.n

    def execute(self):
        all_scores = []
        for genome_index in range(self.n):
            score_objs = self.test_genome(int(genome_index))
            all_scores.append(score_objs)
        return all_scores

    def test_genome(self, genome_index):
        score_objs = []
        i = 0
        genome = self.group.genomes[genome_index]
        while i < self.count_per_g[genome_index]:
            c, rest_g = sample_contig(genome, self.x_st)
            score_obj = score_contig(c, rest_g, genome_index, self)
            score_objs.append(score_obj)
            i+=1
        return score_objs

class ScoreCollection:
    def __init__(self):
        self.genome = None
        self.group = []
        self.other = []
    def __str__(self):
        string = str(self.genome)
        string += "\t"
        for score in self.group:
            string += str(score) + ", "
        string += "\t"
        for other_group in self.other:
            for score in other_group:
                string += str(score) + "; "
            string += ", "
        return string
        

class Score:
    def __init__(self, p_value, species_contig):
        self.p_value = p_value
        self.species_contig = species_contig
        self.genus_contig = ""
        self.species_compare = ""
        self.genus_compare = ""
    def __str__(self):
        return "%f\t%s\t%s\t%s\t%s" % (self.p_value, self.species_contig, self.genus_contig, self.species_compare, self.genus_compare)
