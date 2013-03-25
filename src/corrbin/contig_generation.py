from random import randint
from probin.dna import DNA

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class SampleSetting(object):
    """A class for storing variables related to the contig sampling"""
    def __init__(self, prio, no_contigs, contig_min_length, contig_max_length, debug_mode):
        self.prio = prio
        self.no_contigs = no_contigs
        self.contig_min_length = contig_min_length
        self.contig_max_length = contig_max_length
        self.debug_mode = debug_mode

def sample_contig(genome, x_st, contig_id, start_position=False):
    """ Generates a contig from genome genome

    :genome - DNA object
    :x_st - SampleSetting object
    :contig_id - The unique id given to this contig"""
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
    if start_position:
        contig.start_position = start
    return contig


class SampleGroup(object):
    """Somewhat a duplicate of multinomial.Experiment"""
    def __init__(self,s_st, group,id_gen):
        self.s_st = s_st
        self.group = group
        self.count_per_g = self.count_per_g()
        self.id_gen = id_gen

    def count_per_g(self):
        if self.s_st.prio == "genomes":
            return [round(self.s_st.no_contigs/float(self.n))]*self.n
        elif self.s_st.prio == "groups":
            base = self.s_st.no_contigs/self.n
            count_l = [base]*self.n
            list_index = 0
            while sum(count_l)<self.s_st.no_contigs:
                try:
                    count_l[list_index] += 1
                    list_index += 1
                except:
                    list_index = 0
            return count_l

    @property
    def n(self):
        return len(self.group)

    def generate_group_contigs(self, start_position=False):
        for genome_index in range(self.n):
            genome = self.group.genomes[genome_index]
            genome.contigs = []
            i = 0 
            while i< self.count_per_g[genome_index]:
                c = sample_contig(genome, self.s_st, self.id_gen.id(), start_position=start_position)
                genome.contigs.append(c)
                i+=1
        return None

    def print_group_contigs(self,file_handle, start_position=False):
        recs = []
        for genome in self.group.genomes:
            for contig in genome.contigs:
                if start_position:
                    id_string = id="{0}_{1}_{2}".format(genome.id,contig.contig_id, contig.start_position)
                else:
                    id_string = id="{0}_{1}".format(genome.id,contig.contig_id)
                sequence = SeqRecord(Seq(contig.full_seq), id_string,name="",description="{0}|{1}|{2}".format(genome.family, genome.genus, genome.species))
                recs.append(sequence)
                
        SeqIO.write(recs,file_handle,"fasta")


