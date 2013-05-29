import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from corrbin.score import parse_contig_description
from corrbin.misc import GenomeGroup
from probin.dna import DNA
import pandas as p
import numpy as np
import sys

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
    
def read_FASTA_files(groups, dir_path,dir_structure='tree'):
    cur_dir = os.getcwd()
    if os.path.isfile(dir_path):
        os.chdir(os.path.dirname(dir_path))
        seq_file = os.path.basename(dir_path)
    else:
        os.chdir(dir_path)
    if dir_structure == 'single_fasta_file':
        seq_list = list(SeqIO.parse(seq_file,"fasta"))
        seq_dic = {}
        for seq in seq_list:
            seq_dic[seq.id] = seq
    for group in groups:
        for genome_data in group.genome_data:
            seq_name = genome_data['file_name']
            if dir_structure == 'tree':
                fasta_files = os.listdir(seq_name)
                for fasta_file in fasta_files:
                    genome_file = open(seq_name + '/' + fasta_file)
                    identifier = genome_file.readline()
                    # Only use non-plasmid genomes
                    # Some bacterial genomes contain more than 1 chromosome,  
                    # but assumed not more than 2
                    if identifier.find('plasmid') == -1 and \
                            (identifier.find('complete genome') != -1 or\
                                 identifier.find('chromosome 1') != -1 or\
                                 identifier.find('chromosome I,') != -1):
                        genome_file.close() #Close and reopen the same file
                        genome_file = open(seq_name + '/' + fasta_file)
                        genome_seq = list(SeqIO.parse(genome_file, "fasta"))
                        if len(genome_seq) > 1:
                            sys.stderr.write("Warning! The file " + fasta_file + " in directory " + dir_name + " contained more than one sequence, ignoring all but the first!" + os.linesep)
                        genome = DNA(id = seq_name, seq= str(genome_seq[0].seq))
                        genome.genus = genome_data['genus']
                        genome.species = genome_data['species']
                        genome.family = genome_data['family']
                        group.genomes.append(genome)
                    genome_file.close()
            elif dir_structure == 'single_fasta_file':
                seq = seq_dic[genome_data['file_name']]
                
                genome = DNA(id = seq.id, seq= str(seq.seq))
                genome.genus = genome_data['genus']
                genome.species = genome_data['species']
                genome.family = genome_data['family']
                group.genomes.append(genome)
                
    os.chdir(cur_dir)

def read_FASTA_files_no_groups(meta_genomes, dir_path,dir_structure='tree'):
    cur_dir = os.getcwd()
    if os.path.isfile(dir_path):
        if os.path.dirname(dir_path) != '':
            os.chdir(os.path.dirname(dir_path))
        seq_file = os.path.basename(dir_path)
    else:
        os.chdir(dir_path) 
    if dir_structure == 'single_fasta_file':
        seq_list = list(SeqIO.parse(seq_file,"fasta"))
        seq_dic = {}
        for seq in seq_list:
            seq_dic[seq.id] = seq
    genomes = []
    for genome_data in meta_genomes:
        dir_name = genome_data['file_name']
        if dir_structure == 'tree':
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
        elif dir_structure == 'single_fasta_file':
            seq = seq_dic[genome_data['file_name']]
            
            genome = DNA(id = seq.id, seq= str(seq.seq))
            genome.genus = genome_data['genus']
            genome.species = genome_data['species']
            genome.family = genome_data['family']
            genomes.append(genome)

    os.chdir(cur_dir)
    return genomes

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

def read_contigs_file(open_contigs_file, start_position=False,taxonomy_info=True,filter_dict=None):
    """ Read contigs file generated by generate_contigs script"""
    
    contigs = []
    seqs = list(SeqIO.parse(open_contigs_file, "fasta"))
    for seq in seqs:
        if filter_dict and not filter_dict[seq.id]:
            continue
        if taxonomy_info:
            contig_id_hash = parse_contig_description(seq.description, start_position=start_position)
            contig = DNA(id=contig_id_hash["genome"], seq=str(seq.seq))
            if start_position:
                contig.start_position = contig_id_hash["start_position"]
            contig.family = contig_id_hash["family"]
            contig.genus = contig_id_hash["genus"]
            contig.species = contig_id_hash["species"]
            contig.contig_id = contig_id_hash["contig_id"]
        else:
            contig = DNA(id=seq.id,seq=str(seq.seq))
            contig.contig_id = seq.id
        contigs.append(contig)

    return contigs

def print_contigs_time_series(cs_with_ts,file_handle,sample_headers):
    header_line ="contig_id\tOTU\t" + "\t".join(sample_headers) + os.linesep
    file_handle.write(header_line)
    for c_with_ts in cs_with_ts:
        line = []
        line.append(str(c_with_ts.contig_id))
        line.append(str(c_with_ts.otu))
        line += [str(ts) for ts in c_with_ts.time_series]
        
        out_line = "\t".join(line) + os.linesep

        file_handle.write(out_line)

def read_time_series(ts_file):
    return p.io.parsers.read_table(ts_file,sep='\t')

def read_time_series_file_genomes(genomes,genome_time_series_file,contig_strain_otu_dic,first_data, last_data,tot_num_bases_per_sample):
    from probin.model.coverage import log_coverage
    ts_df = read_time_series(genome_time_series_file)
    
    otu_genome_dic = {}

    
    for genome in genomes:
        strain = genome.species
        otu = contig_strain_otu_dic[strain]
        genome_length= len(genome.full_seq)
        rf = np.array(ts_df[ts_df['# OTU'] == otu].ix[:,first_data:last_data].values)
        if rf.shape[0]==0:
            sys.stderr.write(otu+'\n')
        if rf.shape[0]>=2:
            sys.stderr.write(otu+' More than 2 \n')
        
        genome.mapping_reads = log_coverage.relative_frequency_to_log_coverage(rf,genome_length,tot_num_bases_per_sample/genome_length)
        
