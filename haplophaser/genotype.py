"""
	genotype.py

	Genotype

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

class Genotype(object):

	# Encode pairs of alleles as 0s, 1s, and 2s
	HOMO_REF = 0
	HOMO_ALT = 1
	HETERO = 2

	sequence = []
	tag = None

	def __init__(self, sequence):
		self.sequence = sequence

		for snp in self.sequence:
			if not (snp == Genotype.HOMO_REF or snp == Genotype.HOMO_ALT
				or snp == Genotype.HETERO):
				raise ValueError("Invalid genotype sequence.")

		return


	# Tags the genotype with an identifier
	def setTag(self, tag):
		self.tag = tag


	def __eq__(self, other):
		if isinstance(other, Genotype):
			return self.sequence == other.sequence

		return False

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash(tuple(self.sequence))

	def __len__(self):
		return len(self.sequence)

	def __str__(self):
		s = ""
		for snp in self.sequence:
			s += str(snp)
		return s
