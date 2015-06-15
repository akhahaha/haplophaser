"""
	haplophaser.py

	Haplotype phase runner

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

from genotype import Genotype
from haplotype import Haplotype
from phaseMatch import PhaseMatch
from greedyPhaser import GreedyPhaser
from spockPhaser import SpockPhaser

from itertools import product
from random import randrange
import time


class Haplophaser:
	""" Haplotype phaser """

	# Algorithm arg values
	GREEDY = 0
	SPOCK = 1

	# Output arg values
	DEFAULT = 0		# Prints run statistics only
	SUPPRESS = 1	# Does not output any text
	VERBOSE = 2		# Prints out phase results

	MAX_LENGTH = 15

	N = None
	M = None
	gtypes = None
	htypes = None

	# TODO Create main functionality
	# def main()
	#
	# if __name__ == '__main__':
	#	main

	def __init__(self):
		return

	# Generate a sample of N genotypes with M SNPs each.
	def generate(self, N, M):
		self.N = N
		self.M = M
		self.gtypes = [None] * N
		for i in xrange(self.N):
			self.gtypes[i] = Genotype([randrange(3) for j in xrange(self.M)])
			self.gtypes[i].setTag(i)
		return self.gtypes


	# Generate a NxM each with P unique haplotype pairs
	def generatePopulation(self, N, M, P):
		self.N = N
		self.M = M
		self.gtypes = [None] * N
		P = min(pow(2, M), P) # Max number of unique pairs is 2^M

		# Cross random haplotypes to create P unique pairs
		gtypes = set()
		while len(gtypes) < P:
			h1 = Haplotype([randrange(2) for j in xrange(self.M)])
			h2 = Haplotype([randrange(2) for j in xrange(self.M)])
			gtypes.add(h1.getGenotype(h2, len(gtypes)))
		gtypes = list(gtypes)

		# Select N genotypes out of the generated unique pairs
		for i in xrange(self.N):
			self.gtypes[i] = gtypes[randrange(P)-1]
			self.gtypes[i].index = i

		return self.gtypes


	# Load sample from file.
	def load(self, fname):
		self.fname = fname;
		f = open(fname, "r")

		# Read in file parameters
		fparams = f.readline().split(" ")
		self.N = int(fparams[0])
		self.M = int(fparams[1])

		# Read in genotypes
		self.gtypes = [None] * self.N
		for i in xrange(self.N):
			genotype = [None] * self.M
			for j in xrange(self.M):
				snp = f.read(1)
				if (snp.isdigit()):
					genotype[j] = int(snp)
			f.read(1) # Read newline
			self.gtypes[i] = Genotype(genotype, i)

		f.close()
		return


	# Save sample to file.
	# TODO def toFile(self, fname):

	# Print current sample information.
	def __str__(self):
		s = "N: " + str(self.N) + " M: " + str(self.M) + "\n"
		for gtype in self.gtypes:
			s += str(gtype) + "\n"
		return s

	# Run current sample with specified algorithm
	def run(self, argAlg, argOutput):
		results = None

		# Set algorithm
		algName = ""
		if (argAlg == self.GREEDY):
			algName = "GREEDY"
		elif (argAlg == self.SPOCK):
			algName = "SPOCK"
		else:
			raise ValueError("Invalid argument: algorithm")

		# Run selected algorithm
		start = time.time()
		if (argAlg == self.GREEDY):
			results = self.runGreedy()
		elif (argAlg == self.SPOCK):
			results = self.runSpock()
		else:
			raise ValueError("Invalid argument: algorithm")
		runtime = time.time() - start

		# Process results
		# TODO Check unique pairs with sample pairs.
		count = len(results)
		uniquePairs = self.countUniquePairs(results)
		if (argOutput == self.VERBOSE):
			self.printPhaseMatches(results)

		# Print run statistics
		if (argOutput != self.SUPPRESS):
			print str(self.N) + "x" + str(self.M) + "\t" + algName + \
				"\t Solution Size: " + str(count) + "\t" + str(uniquePairs) + \
				"\t Time: " + str(runtime)

		return results


	def countUniquePairs(self, phaseMatches):
		uniquePairs = 0
		for pm in phaseMatches:
			uniquePairs += pm.getUniquePairCount()
		return uniquePairs


	def printPhaseMatches(self, phaseMatches):
		for pm in phaseMatches:
			print pm
		return


	# Run greedy phasing algorithm
	def runGreedy(self):
		phaser = GreedyPhaser()
		results = phaser.phase(self.M, self.gtypes)
		return results


	# Run Spock's phasing algorithm
	def runSpock(self):
		phaser = SpockPhaser()

		if (self.M <= self.MAX_LENGTH):
			return phaser.phase(self.M, self.gtypes)

		# Calculate number of windows
		numWindows =  self.M / self.MAX_LENGTH
		if self.M % self.MAX_LENGTH > 0:
			numWindows += 1

		N = len(self.gtypes)
		resultHtypes = [Haplotype([]) for _ in xrange(N)]
		hDict = {}

		for i in xrange(numWindows):
			# Create split matrix
			winG = [None] * N
			for j in xrange(N):
				snip = self.gtypes[j].sequence[i*self.MAX_LENGTH :(i+1)*self.MAX_LENGTH]
				winG[j] = Genotype(snip, j)

			# Process window
			winResults = phaser.phase(len(winG[0]), winG)
			# Process window results
			for pm in winResults:
				for g in pm.gtypes:
					resultHtypes[g.index].sequence += pm.htype.sequence

		# Compile haplotype dictionary
		for i in xrange(N):
			hstr = str(resultHtypes[i])
			if not (hstr in hDict):
				h = PhaseMatch(resultHtypes[i])
			else:
				h = hDict[hstr]
			h.addGenotype(self.gtypes[i])
			hDict[hstr] = h

		return list(hDict.values())
