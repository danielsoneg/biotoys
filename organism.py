from genome import Genome

class Organism:
    @classmethod
    def as_offspring(cls, parent):
        child_no = len(parent.children)
        return cls(parent.genome.as_offspring(parent.genome, child_no), parent, child_no)

    def __init__(self, name, genome, parent=None, child_no=0):
        self.genome = genome
        self.parent = parent
        self.child_no = child_no
        self.children = []

    @classmethod
    def new(cls, genome_length):
        return cls(Genome.generate(genome_length))

    def __lineage(self):
        if self.parent:
            yield from self.parent.__lineage
        yield self.parent

    @property
    def siblings(self):
        return self.parent.children

    @property
    def lineage(self):
        return " <- ".join(self.__lineage)
    
    #def __repr__(self):
    
    def reproduce(self, count):
        for i in range(count):
            self.children.append(Organism.as_offspring(self))
