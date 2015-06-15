"""
	haplotype.py

	Haplotype

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

from genotype import Genotype


class Haplotype(object):

	# Encode alleles as 0s and 1s
	REF = 0
	ALT = 1

	def __init__(self, sequence):
		self.sequence = sequence

		for snp in self.sequence:
			if not (snp == Haplotype.REF or snp == Haplotype.ALT):
				raise ValueError("Invalid haplotype sequence.")

		return

	# Determines whether the haplotype explains a given genotype.
	def match(self, genotype):
		gtype = genotype.sequence
		if (len(gtype) != len(self.sequence)):
			raise ValueError("Genotype/Haplotype length mismatch")

		for i in xrange(len(gtype)):
			if not (self.sequence[i] == Haplotype.REF
				or self.sequence[i] == Haplotype.ALT):
				raise ValueError("Invalid haplotype sequence.")
			if not (gtype[i] == Genotype.HETERO
				or gtype[i] == self.sequence[i]):
				return False
		return True


	# Generates a complement for the genotype/haplotype.
	def getComplement(self, genotype):
		gtype = genotype.sequence
		if (len(gtype) != len(self.sequence)):
			raise ValueError("Genotype/Haplotype length mismatch")

		complement = [None] * len(gtype)
		for i in xrange(len(self.sequence)):
			if not (self.sequence[i] == Haplotype.REF
				or self.sequence[i] == Haplotype.ALT):
				raise ValueError("Invalid haplotype sequence.")

			if (gtype[i] == Genotype.HETERO):
				# Ambiguous SNP
				complement[i] = Haplotype.REF if (self.sequence[i] == \
					Haplotype.ALT) else Haplotype.ALT
			elif (gtype[i] == Genotype.HOMO_REF
				or gtype[i] == Genotype.HOMO_ALT):
				# Unambiguous SNP
				complement[i] = gtype[i]
			else:
				raise ValueError("Invalid genotype sequence.")

		return Haplotype(complement)


	# Returns the genotype when paired with another haplotype
	def getGenotype(self, haplotype, index):
		htype = haplotype.sequence
		if (len(htype) != len(self.sequence)):
			raise ValueError("Haplotype/Haplotype length mismatch")

		gtype = [None] * len(htype)
		for i in xrange(len(self.sequence)):
			if (self.sequence[i] == htype[i]):
				if not (self.sequence[i] == Haplotype.REF
				or self.sequence[i] == Haplotype.ALT):
					raise ValueError("Invalid haplotype sequence.")

				gtype[i] = self.sequence[i]
			else:
				gtype[i] = Genotype.HETERO

		return Genotype(gtype, index)


	def __eq__(self, other):
		if isinstance(other, Haplotype):
			return self.sequence == other.sequence

		return False

	def __ne__(self, other):
		return not self.__eq__(other)

	# Returns the length of the haplotype.
	def __len__(self):
		return len(self.sequence)

	# Returns the haplotype string.
	def __str__(self):
		s = ""
		for snp in self.sequence:
			s += str(snp)
		return s
