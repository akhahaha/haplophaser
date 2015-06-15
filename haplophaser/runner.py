"""
	runner.py

	Helper functions to run Haplophaser.

	Author(s):
		Alan Kha		akhahaha@gmail.com

"""

import genotype
import haplotype
import phaseMatch
import greedyPhaser
import spockPhaser
import haplophaser
from haplophaser import Haplophaser
import runner

FIXED = 0
SIZE = 1
LENGTH = 2

hp = Haplophaser()

def refresh():
	reload(genotype)
	reload(haplotype)
	reload(phaseMatch)
	reload(greedyPhaser)
	reload(spockPhaser)
	reload(haplophaser)
	reload(runner)
	print "All dependent modules reloaded"
	return

# Runs an algorithm against a generated NxM sample
def run(N, M, flag, alg):
	if (flag == FIXED):
		hp.generate(N, M)
		hp.run(alg, Haplophaser.DEFAULT)
	elif (flag == SIZE):
		for i in xrange(N):
			hp.generate(i+1, M)
			hp.run(alg, Haplophaser.DEFAULT)
	elif (flag == LENGTH):
		for i in xrange(M):
			hp.generate(N, i+1)
			hp.run(alg, Haplophaser.DEFAULT)
	return

def runCompare(N, M, flag):
	if (flag == FIXED):
		hp.generate(N, M)
		hp.run(0, Haplophaser.DEFAULT)
		hp.run(1, Haplophaser.DEFAULT)
	elif (flag == SIZE):
		for i in xrange(N):
			hp.generate(i+1, M)
			hp.run(0, Haplophaser.DEFAULT)
			hp.run(1, Haplophaser.DEFAULT)
	elif (flag == LENGTH):
		for i in xrange(M):
			hp.generate(N, i+1)
			hp.run(0, Haplophaser.DEFAULT)
			hp.run(1, Haplophaser.DEFAULT)
	return

def runT(N, M, alg):
	hp.generate(N, M)
	print str(hp)
	hp.run(alg, Haplophaser.VERBOSE)
	return

def runP(N, M, P, alg):
	hp.generatePopulation(N, M, P)
	hp.run(alg, Haplophaser.DEFAULT)
	return

def runTP(N, M, P, alg):
	hp.generatePopulation(N, M, P)
	hp.run(alg, Haplophaser.VERBOSE)
	return

def runPSet(N, M, Pset, alg):
	for P in Pset:
		hp.generatePopulation(N, M, P)
		hp.run(alg, Haplophaser.DEFAULT)
	return

def runTF(F, alg):
	hp.load(F)
	hp.run(alg, Haplophaser.VERBOSE)
	return

set1 = [1, 10, 25, 50, 100, 200, 500, 750, 1000, 1500, 2000]
def runSet(Nset, Mset, alg):
	for N in Nset:
		for M in Mset:
			hp.generate(N, M)
			hp.run(alg, Haplophaser.DEFAULT)
	return
