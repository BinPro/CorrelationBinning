from probin.dna import DNA
from random import randint
from corrbin.misc import GenomeGroup
import os
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

def sample_contig(genome, x_st, contig_id):
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
    
    return contig

def read_parsed_taxonomy_file(open_name_file):
    groups = []
    # Read the file with all names, divide them into groups
    for line in open_name_file:
        if line[0:12] == 'family_name:':
            family = line.split('\t')[1].strip()
        elif line[0:11] == 'genus_name:':
            genus = line.split('\t')[1].strip()
            new_group = GenomeGroup(genus)
            new_group.family = family
            groups.append(new_group)
        elif line[0:6] == 'entry:':
            genome_name = line.split('\t')[2].strip()
            genome_species = line.split('\t')[1].strip()
            meta_genome = {'id': genome_name,
                           'species': genome_species,
                           'genus': genus,
                           'family': family,
                           'file_name': genome_name
                          }
            groups[-1].genome_data.append(meta_genome)
    return groups

def read_FASTA_files(groups, dir_path):
    cur_dir = os.getcwd()
    os.chdir(dir_path)
    for group in groups:
        for genome_data in group.genome_data:
            dir_name = genome_data['file_name']
            fasta_files = os.listdir(dir_name)
            for fasta_file in fasta_files:
                genome_file = open(dir_name + '/' + fasta_file)
                identifier = genome_file.readline()
                # Only use non-plasmid genomes
                # Some bacterial genomes contain more than 1 chromosonme,  
                # but assumed not more than 2
                if identifier.find('plasmid') == -1 and identifier.find('chromosome 2') == -1:
                    genome_file.close() #Close and reopen the same file
                    genome_file = open(dir_name + '/' + fasta_file)
                    genome_seq = list(SeqIO.parse(genome_file, "fasta"))
                    if len(genome_seq) > 1:
                        sys.stderr.write("Warning! The file " + fasta_file + " in directory " + dir_name + " contained more than one sequence, ignoring all but the first!" + os.linesep)
                    genome = DNA(id = dir_name, seq= str(genome_seq[0].seq))
                    genome.genus = genome_data['genus']
                    genome.species = genome_data['species']
                    genome.family = genome_data['family']
                    group.genomes.append(genome)
                genome_file.close()

    os.chdir(cur_dir)

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

    def generate_group_contigs(self):
        for genome_index in range(self.n):
            genome = self.group.genomes[genome_index]
            genome.contigs = []
            i = 0 
            while i< self.count_per_g[genome_index]:
                c = sample_contig(genome, self.s_st, self.id_gen.id())
                genome.contigs.append(c)
                i+=1
        return None

    def print_group_contigs(self,file_handle):
        recs = []
        for genome in self.group.genomes:
            for contig in genome.contigs:
                sequence = SeqRecord(Seq(contig.full_seq), id="{0}_{1}".format(genome.id,contig.contig_id), name="",description="{0}|{1}".format(genome.family, genome.genus))
                recs.append(sequence)
                
        SeqIO.write(recs,file_handle,"fasta")
