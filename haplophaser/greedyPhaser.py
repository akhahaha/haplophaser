"""
	greedyPhaser.py

	Greedy haplotype phaser

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

from phaser import Phaser

from genotype import Genotype
from haplotype import Haplotype
from phaseMatch import PhaseMatch

from itertools import product


class GreedyPhaser(Phaser):
	"""
		Greedy haplotype phaser

		The greedy phasing algorithm generates all 2^M possible haplotypes and
		selects the haplotype that explains the most unresolved genotypes and
		repeats until all genotypes have been resolved.

	"""

	def __init__(self):
		return


	def phase(self, M, gtypes):
		N = len(gtypes)
		if (N == 0):
			return None

		# Generate all possible haplotypes of length M
		resolved = 0
		htypes = self.generateHaplotypes(M)
		results = []

		while resolved < N:
			# Find the haplotype that explains the most unresolved genotypes
			best = None
			bestCount = 0
			bestMatches = [False] * len(gtypes)
			possible = []

			for i in xrange(len(htypes)):
				count = 0
				matches = [False] * len(gtypes)

				for j in xrange(len(gtypes)):
					if (M != len(gtypes[j])):
						raise ValueError("Incorrect genotype length.")
					if htypes[i].match(gtypes[j]):
						count += 1
						matches[j] = True

				if count > bestCount:
					if not (best == None):
						possible.append(best)
					best = htypes[i]
					bestCount = count
					bestMatches = matches
				elif count != 0:
					possible.append(htypes[i])
				# Discard haplotypes with no matches

			# Resolve best matches
			pm = PhaseMatch(best)
			unresolved = []
			for i in xrange(len(gtypes)):
				if bestMatches[i]:
					pm.addGenotype(gtypes[i])
					resolved += 1
				else:
					unresolved.append(gtypes[i])
			results.append(pm)
			gtypes = unresolved

			# Use reduced set of possible haplotypes (-no matches, -last best)
			htypes = possible

			if len(possible) == 0:
				break

		return results


	def generateHaplotypes(self, M):
		""" Generates all possible haplotypes of length M """

		htypes = []
		for hsequence in product([0, 1], repeat = M):
			htypes.append(Haplotype(hsequence))

		return htypes
