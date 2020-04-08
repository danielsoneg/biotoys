import random
import functools
import string

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VALID = ("A", "T", "G", "C")
MASK = "N"
r_mutate  = 0.0
r_crop    = 0.0
r_degrade = 0.0
r_snp     = 0.0

def set_params(mutate, crop, degrade, snp):
    global r_mutate, r_crop, r_degrade, r_snp
    r_mutate, r_crop, r_degrade, r_snp = mutate, crop, degrade, snp

class Sequence:
    __base_seq = ""
    name = ""

    @classmethod
    def generate(cls, length, name=None):
        name = name or "".join(random.choices(string.ascii_uppercase, k=3))
        return cls(name, generate(length))
    
    @classmethod
    def as_offspring(cls, seq, index=1):
        return cls(seq.name + ".%s" % index, seq.mutate())

    def __init__(self, name, seq):
        self.name = name
        self.base_seq = seq
    
    def __repr__(self):
        return "Sequence %s: %s" % (self.name, self.seq)

    def mutate(self):
        return mutate(r_mutate, self.base_seq, r_snp=r_snp)

    @property
    def base_seq(self):
        return self.__base_seq

    @base_seq.setter
    def base_seq(self, seq):
        if hasattr(self, "__seq"):
            del(self.__seq)
        self.__base_seq = seq

    @property
    def seq(self):
        if not hasattr(self, "__seq"):
            self.__seq = crop(r_crop,
                         degrade(r_degrade,
                         self.__base_seq))
        return self.__seq
    
    def reproduce(self, offspring):
        return [Sequence.as_offspring(self, index=s) for s in range(offspring)]
    
    @property
    def seq_record(self):
        return SeqRecord(Seq(self.seq), id=self.name, description=self.description)

    @property
    def generation(self): return self.name.count(".")
    @property
    def childnumber(self): return self.name.rpartition(".")[2] 
    @property
    def lineage(self): return self.name.partition(".")[0]
    @property
    def description(self):
        return "Child %s of generation %s of lineage %s" % (
                self.childnumber, self.generation, self.lineage)

def generate(length):
    return "".join(random_nucleotide() for i in range(length))

def mutate(rate, seq, r_snp=0):
    seq = (
        base if random.random() >= rate else mutate_nucleotide(base)
        for base in seq
    )
    if r_snp:
        seq = snp(r_snp, seq)
    return "".join(seq)

def snp(count, seq):
    seq = list(seq)
    for idx in random.sample(range(len(seq)), k=count):
        existing = seq[idx]
        while seq[idx] == existing: # No for real mutate
            seq[idx] = random_nucleotide()
    return list(seq)

def random_nucleotide():
    return random.choice(VALID)

def mutate_nucleotide(base):
    return random.choice((
        "",
        random_nucleotide(),
        "".join(random_nucleotide() for r in range(random.randint(0,3)))
    ))

def crop(rate, seq):
    mu = len(seq) * rate
    sigma = mu * .2
    drop = int(random.gauss(mu, sigma))
    beginning = int(drop * random.random())
    end = drop - beginning
    return seq[beginning:-end or None]

def degrade(rate, seq):
    return "".join(base if random.random() >= rate else "N" for base in seq)

def generate_tree(seq_length, generations, rate, name=None):
    tree = list(iter_tree(Sequence.generate(seq_length, name=name), rate, 0, generations))
    print("Root: %s" % tree[0].name)
    return tree

def iter_tree(seq, per_gen, current, max_gen):
    yield seq
    if current < max_gen:
        for child in seq.reproduce(per_gen):
            yield from iter_tree(child, per_gen, current + 1, max_gen)

def output_tree(tree, dest):
    with open(dest, "w") as fasta:
        SeqIO.write((s.seq_record for s in tree), fasta, "fasta")