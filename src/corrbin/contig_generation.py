from probin.dna import DNA
from random import randint


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

def read_FASTA_files(groups):
    pass
