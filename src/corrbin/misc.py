class Uniq_id(object):
    def __init__(self, start_val=0, incr=1):
        self.start = start_val
        self.curr_val = start_val
        self.incr = incr

    def id(self):
        val = self.curr_val
        self.curr_val +=self.incr
        return val

class GenomeGroup(object):
    """A wrapper for the genome groups"""
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

def all_but_index(l,i):
    return l[0:i] + l[i+1:]
