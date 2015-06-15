"""
	spockPhaser.py

	Spock haplotype phaser

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

from phaser import Phaser

from genotype import Genotype
from haplotype import Haplotype
from phaseMatch import PhaseMatch

from copy import deepcopy
from itertools import product


class SpockPhaser(Phaser):
	"""
		Spock haplotype phaser

		The Spock phasing algorithm uses a depth-first search branching on the
		SNP with the most information (lowest heterozygous frequency),
		and using the more common homozygous value at that position.
	"""

	def __init__(self):
		return


	def phase(self, M, gtypes):
		bestCount = len(gtypes) # Max solution size is N
		bestResults = []

		# Run initial deduction
		dStack = [self.deduce(M, deepcopy(gtypes), [], 0, bestCount)]

		# Depth-first search to find the optimal solution
		while dStack:
			d = dStack.pop()

			if d.results and not d.unresolved:
				if d.resultCount <= bestCount:
					bestCount = d.resultCount
					bestResults = d.results
			elif d.unresolved:
				# Continue tree
				urREF = d.unresolved
				urALT = d.unresolved

				for i in xrange(len(d.unresolved)):
					# Branch factor: highest frequency heterogeneous SNP
					if (d.unresolved[i].sequence[d.nextBranch] == Genotype.HETERO):
						urREF[i].sequence[d.nextBranch] = Genotype.HOMO_REF
						urALT[i].sequence[d.nextBranch] = Genotype.HOMO_ALT

				dREF = self.deduce(M, urREF, d.results, d.resultCount, bestCount)
				dALT = self.deduce(M, urALT, d.results, d.resultCount, bestCount)

				# Determine branch order
				if (d.nextValue == Genotype.HOMO_REF):
					dStack.append(dALT)
					dStack.append(dREF)	# Search REF first
				else:
					dStack.append(dREF)
					dStack.append(dALT)	# Search ALT first

		# Restore original genotype sequences
		for i in xrange(len(bestResults)):
			p = PhaseMatch(bestResults[i].htype)
			for gtype in bestResults[i].gtypes:
				p.addGenotype(gtypes[gtype.tag])
			bestResults[i] = p

		return bestResults


	def deduce(self, M, gtypes, prevResults, prevCount, pruneCount):
		"""
			Eliminate resolved, and determine the best SNP and value to perform
			the next branch with.
		"""

		# TODO Use global frequencies instead of per deduction

		N = len(gtypes)

		unresolved = []
		resolved = []

		diffMap = [0] * M

		freq = [0] * M
		mode = None
		modeVal = 0

		results = []

		for i in xrange(N):
			gtype = gtypes[i]

			if (M != len(gtype)):
				raise ValueError("Incorrect genotype length.")
			count = 0

			rCandidates = results
			for j in xrange(M):
				if (gtype.sequence[j] == Genotype.HETERO):
					freq[j] += 1
					count += 1
					if (freq[j] > modeVal):
						mode = j
						modeVal = freq[j]
				# Eliminate resolved candidates on SNP contradiction
				else:
					if (gtype.sequence[j] == Genotype.HOMO_REF):
						diffMap[j] += 1
					elif (gtype.sequence[j] == Genotype.HOMO_ALT):
						diffMap[j] -= 1

				# Check for resolved haplotypes
				rCandidates = [rc for rc in rCandidates if not (
					(rc.htype.sequence[j] != gtype.sequence[j]))]

			# Already resolved/deduced haplotype
			if (rCandidates):
				# Add to earliest result
				results[results.index(rCandidates[0])].addGenotype(gtype)
				resolved.append(gtype)
			# Newly resolved haplotype
			elif (count == 0):
				pm = PhaseMatch(Haplotype(gtype.sequence))
				pm.addGenotype(gtype)
				resolved.append(gtype)
				results.append(pm)

				# Prune if necessary
				if prevCount + len(resolved) > pruneCount:
					return Deduction(M, [], 0, gtypes, -1, -1)
			# Unresolved haplotype
			else:
				unresolved.append(gtype)

		# Calculate final heterogeneous SNP frequency
		for i in xrange(M):
			for gtype in resolved:
				if (gtype.sequence[j] == Genotype.HETERO):
					freq[j] -= 1
					if (j == mode):
						modeVal -= 1
					elif (freq[j] > modeVal):
						mode = j
						modeVal = freq[j]
				elif (gtype.sequence[j] == Genotype.HOMO_REF):
					diffMap[j] -= 1
				elif (gtype.sequence[j] == Genotype.HOMO_ALT):
					diffMap[j] += 1

		# Difference map heuristic
		# Eliminate positions with no heterogeneous SNPs
		for i in xrange(M):
			if freq[i] == 0:
				diffMap[i] = 0

		modeRefVal = max(diffMap)
		modeAltVal = min(diffMap)
		if (modeRefVal > abs(modeAltVal)):
			modeDiffVal = modeRefVal
		else:
			modeDiffVal = modeAltVal
		modeDiffPos = diffMap.index(modeDiffVal)

		# Ensure chosen value has no change
		# TODO Choose next highest value
		if (modeDiffPos != 0):
			nextValue = modeDiffPos
			nextGuess = 1 if modeDiffVal > 0 else 0
		else:
			nextValue = mode
			nextGuess = Genotype.HOMO_REF # Default to reference allele

		return Deduction(M, results+prevResults, prevCount + len(resolved),
			unresolved, nextValue, nextGuess)



class Deduction:
	""" Deduce result object """

	def __init__(self, M, results, resultCount, unresolved, nextBranch, nextValue):
		self.M = M
		self.results = results
		self.resultCount = resultCount
		self.unresolved = unresolved
		self.nextBranch = nextBranch
		self.nextValue = nextValue
		return
