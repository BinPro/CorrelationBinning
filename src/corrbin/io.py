import os
from Bio import SeqIO

from corrbin.misc import GenomeGroup
from probin.dna import DNA


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

def genome_info_from_parsed_taxonomy_file(open_name_file):
    genomes = []
    for line in open_name_file:
        if line[0:12] == 'family_name:':
            family = line.split('\t')[1].strip()
        elif line[0:11] == 'genus_name:':
            genus = line.split('\t')[1].strip()
        elif line[0:6] == 'entry:':
            genome_name = line.split('\t')[2].strip()
            genome_species = line.split('\t')[1].strip()
            meta_genome = {'id': genome_name,
                           'species': genome_species,
                           'genus': genus,
                           'family': family,
                           'file_name': genome_name
                          }
            genomes.append(meta_genome)
    return genomes
    
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
                # Some bacterial genomes contain more than 1 chromosome,  
                # but assumed not more than 2
                if identifier.find('plasmid') == -1 and \
                        (identifier.find('complete genome') != -1 or\
                             identifier.find('chromosome 1') != -1 or\
                             identifier.find('chromosome I,') != -1):
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

def read_FASTA_files_no_groups(meta_genomes, dir_path):
    cur_dir = os.getcwd()
    os.chdir(dir_path)
    genomes = []
    for genome_data in meta_genomes:
        dir_name = genome_data['file_name']
        fasta_files = os.listdir(dir_name)
        for fasta_file in fasta_files:
            genome_file = open(dir_name + '/' + fasta_file)
            identifier = genome_file.readline()
            # Only use non-plasmid genomes
            # Some bacterial genomes contain more than 1 chromosome,
            # but assumed not more than 2
            if identifier.find('plasmid') == -1 and \
                    (identifier.find('complete genome') != -1 or\
                         identifier.find('chromosome 1') != -1):
                # Close and reopen the same file
                genome_file.close()
                genome_file = open(dir_name + '/' + fasta_file)
                genome_seq = list(SeqIO.parse(genome_file, "fasta"))
                if len(genome_seq) > 1:
                    sys.stderr.write("Warning! The file " + fasta_file + " in directory " + dir_name + " contained more than one sequence, ignoring all but the first!" + os.linesep)
                genome = DNA(id = dir_name, seq= str(genome_seq[0].seq))
                genome.genus = genome_data['genus']
                genome.species = genome_data['species']
                genome.family = genome_data['family']
                genomes.append(genome)
            genome_file.close()
    os.chdir(cur_dir)
    return genomes


def read_genome_parts(file_handle):
    genomes = {}
    for part in SeqIO.parse(file_handle, "fasta"):
        genome = part.id.split("_")

def print_parts(parts,file_handle, ids,genome):
    recs =[]
    for part_number, part in enumerate(parts):
        uniq_id = ids.id
        sequence = SeqRecord(\
            Seq(part.full_seq), 
            id="{0}_{1}_{2}".format(genome.id,
                                    uniq_id(),
                                    part.start_position),
            name="",
            description="{0}|{1}|{2}".format(genome.family, 
                                             genome.genus, 
                                             genome.species, 
                                             part.start_position))
        recs.append(sequence)
    SeqIO.write(recs,file_handle,"fasta")
