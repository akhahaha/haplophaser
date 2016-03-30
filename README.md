Haplophaser
===================
	UCLA CS CM124 Final Project
	Spring 2015

	Alan Kha        904030522	akhahaha@gmail.com
-------------------------------------------------------------------------------
Overview
---------------
Minimum parsimony haplotype phaser.

Features
---------------
 - Population generation
 - Import population data from file
 - Multiple phasing algorithms
   1. Greedy algorithm
   2. Spock's algorithm (DFS branching on most-information SNPs)

Data Format Example
---------------
    7 6
    222222
    000001
    022002
    222222
    000001
    021002
    000001

Installation
---------------
1. Use helper functions in [runner.py](./haplophaser/runner.py) or call
[Haplophaser](./haplophaser/haplophaser.py) directly.
