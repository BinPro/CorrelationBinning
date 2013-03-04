from Bio import SeqIO, SeqRecord, Seq

def get_sequences(file_path):
    seqs = SeqIO.parse(file_path,"fasta")
    return list(seqs)

def get_contigs(seqs,length):
    d = {}
    for seq in seqs:
        d[seq.id] = [seq.seq[i:i+length] for i in xrange(0,len(seq.seq),length)]
        if len(d[seq.id][-1]) != length:
            d[seq.id].pop()
    return d

def write_sequences(seqs, file_path):
    recs = []
    for genome in sorted(seqs.keys()):
        contigs = seqs[genome]
        for i in xrange(len(contigs)):
            recs.append(SeqRecord.SeqRecord(contigs[i],id="{0}_{1}".format(genome,i),name="",description=""))
    SeqIO.write(recs,file_path,"fasta")

if __name__=="__main__":
    import sys
    seqs = get_sequences(sys.argv[1])
    contigs = get_contigs(seqs,sys.argv[2])
    write_sequences(seqs, sys.argv[3])