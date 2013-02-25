from copy import copy
from random import randint

from probin.dna import DNA
from probin.model.composition import multinomial as mn

from corrbin.score import ScoreCollection, Score

def sample_contig(genome, x_st, contig_id):
    min_length = x_st.contig_min_length
    max_length = x_st.contig_max_length
    l = randint(min_length, max_length)
    gen_l = len(genome.full_seq)
    if x_st.debug_mode:
        start = 0
    else:
        start = randint(0, (gen_l-l))
    end = start+l
    contig = DNA(id = genome.id + " contig", seq = genome.full_seq[start:end])
    contig.contig_id = contig_id
    
    rest_par = par_subtraction(copy(genome.signature),contig) 
    return contig, rest_par

def par_subtraction(signature, contig):
    signature.subtract(contig.signature)
    return signature_to_parameters(signature)

def signature_to_parameters(signature):
    pa = {}
    n = sum(signature.values())
    for i in xrange(DNA.kmer_hash_count):
        pa[i] = signature[i]/float(n)
    return pa

def score_contig(contig, par, genome, genome_index, test):
    score_obj = ScoreCollection()
    score_gen = Score(mn.log_probability(contig.signature, par), genome, genome, contig.contig_id)
    score_obj.genome = score_gen
    group_genomes = test.group.all_genomes_but_index(genome_index)
    for gen in group_genomes:
        score_obj.group.append(Score(mn.log_probability(contig.signature, gen.par()),genome, gen, contig.contig_id))
    for group in test.rest_groups:
        outside_group = []
        for gen in group.genomes:
            score = Score(mn.log_probability(contig.signature,gen.par()), genome, gen, contig.contig_id)
            outside_group.append(score)
        score_obj.other.append(outside_group)
    return score_obj

# Quite ugly way of adding a new function to DNA class
# Will be added as soon as any other function from this
# module is loaded...
def par(self):
    try:
        return self.parameters
    except:
        self.parameters = signature_to_parameters(self.signature)
        return self.parameters

DNA.par = par


class GenomeGroup(object):
    """A wrapper for the genome groups used in pairwise testing"""
    def __init__(self, name):
        self.genomes = []
        self.name = name
        self.genome_data = []
    def add_genome(self, genome):
        self.genomes.append(genome)
    def __len__(self):
        return len(self.genomes)
    def all_genomes_but_index(self,i):
        return self.genomes[0:i] + self.genomes[i+1:]

class ExperimentSetting(object):
    """A class for storing variables related to all experiments"""
    def __init__(self, prio, mode, no_contigs, contig_min_length, contig_max_length, debug_mode):
        self.prio = prio
        self.mode = mode
        self.no_contigs = no_contigs
        self.contig_min_length = contig_min_length
        self.contig_max_length = contig_max_length
        self.debug_mode = debug_mode

class Test(object):
    def __init__(self,x_st, group, rest_groups, id_gen):
        self.n = len(group)
        self.count_per_g = self.count_per_g(x_st.no_contigs, self.n)
        self.group = group
        self.rest_groups = rest_groups
        self.x_st = x_st
        self.id_gen = id_gen

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
            c, new_par= sample_contig(genome, self.x_st, self.id_gen.id())
            score_obj = score_contig(c, new_par, genome, genome_index, self)
            score_objs.append(score_obj)
            i+=1
        return score_objs
