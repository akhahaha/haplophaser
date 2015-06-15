"""
	phaseMatch.py

	PhaseMatch

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

from copy import deepcopy


class PhaseMatch:
	""" Stores a phased haplotype with its matched genotypes """

	def __init__(self, htype):
		self.htype = htype
		self.gtypes = []

		return

	def addGenotype(self, gtype):
		if (len(self.htype) != len(gtype)):
			raise ValueError("Genotype/haplotype length mismatch.")
		self.gtypes.append(gtype)

		return

	# Counts the number of unique haplotype pairs for genotype sets
	def getUniquePairCount(self):
		# Generate complements
		complements = [self.htype.getComplement(gtype) for gtype in self.gtypes]
		count = 0

		for i in xrange(len(complements)):
			# Only check if unmatched
			if complements[i] != None:
				count += 1
				# Check complement against all other unmatched genotypes
				for j in xrange(i+1, len(complements)):
					if complements[j] != None and complements[i] == complements[j]:
						complements[j] = None

		return count

	def __str__(self):
		s = ""
		for gtype in self.gtypes:
			s += str(gtype) + "\n"
			s += "\t" + str(self.htype) + "\n"
			s += "\t" + str(self.htype.getComplement(gtype)) + "\n"
		return s
