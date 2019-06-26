#!/users/a2e/ldeolive/myEnv/Env1/bin/python

import argparse
#help for users 
parser = argparse.ArgumentParser()
parser.add_argument("N", type = int, help="Nombre de fragments demandés")
parser.add_argument("fichier_chromosome", help ="fichier contenant 3 colonnes : le nom du chromosome, la position start et la position end.")
parser.add_argument("fichier_préférence_nucléotidique", help = "fichier contenant 2 colonnes : le nucléotide A,T,C ou G et la fréquence associée au nucléotide.")
parser.add_argument("fichier_sortie", help = "Nom du fichier créé contenant 3 colonnes : le nom du chromosome, la position start du fragment et la position end.")
args = parser.parse_args()

import pandas as pd
import numpy as np
import os, sys 
from joblib import Parallel , delayed
import gc
import fileinput, re
import random
import pybedtools
from scipy.stats import johnsonsb
import csv
from tqdm import tqdm

#########################################################
##================= Fonctions =========================##
#########################################################

#Reading file : probability of each nucleotide
def ProbNucl (ProbNuclFile) : 
	listeNucl = []
	listeProb = []
	for line in fileinput.input(ProbNuclFile): 
		liste = line.split("\t")
		Nucl = liste[0]
		Prob = liste[1]
		listeNucl.append(Nucl)
		listeProb.append(Prob)
	return listeNucl, listeProb

def Lenchr(LenchrFile) : 
	listeChromPos = []
	for line in fileinput.input(LenchrFile) : 
		liste = line.split("\n")
		listeChromPos.append(liste[0])
	return listeChromPos

def createRangePybedtools(element,chrom):
	''' This function creates a list of IDs (Chr.start.end) from list end ends. It is called when sampling second cuts around first cuts..'''
	return "{} {} {}".format(chrom, str(element-1),str(element))

def createRange(element,chrom):
	''' This function creates a list of IDs (Chr.start.end) from list end ends. It is called when sampling second cuts around first cuts..'''
	return "{}	{}".format(chrom,str(element))

def SampleElem (array, probs):
	'''Sample one element from an array using probabilities in probs. 
	Here the array is each nucs (for sampling a nucleotide) or list of positions'''
	return np.random.choice(array, 1, p=probs)

#Reading file : length of chromosome
def Lenchr(LenchrFile) : 
	listeChromPos = []
	for line in fileinput.input(LenchrFile) : 
		liste = line.split("\n")
		listeChromPos.append(liste[0])
	return listeChromPos

#tester si une variable est un nombre
def isfloat(str):
	try:
		float(str)
	except ValueError:
		return False
	return True



def getsecondCut (LenEnd, firstCut,fragLens):
	#split first cut
	#print(firstCut)
	if isfloat(firstCut) == True : 
		erange=int(firstCut)
		ends=erange+fragLens # to save time, restrict range to those positions that have a chance to be sampled (freq>0) 
		if ends > LenEnd : 
			ends = LenEnd
		rangeIDs= ends
	else : 
		rangeIDs = "NA"
	
	return rangeIDs

def JohnsonSBdistrib (a , b) : 
	r = johnsonsb.rvs(a, b,  loc=24, scale=1900, size=100000)
	for Lenreads in r : 
		if 30<= Lenreads <= 2000: 
			Lenreads = r
			return r

def fragLen (r) : 
		LenFrag = random.choice(r)
		return LenFrag

def writeFinalCuts (firstCut, secondCut, outFile):
	outputfile = open('{}'.format(outFile),'a')
	writer = csv.writer(outputfile, delimiter = '\t', lineterminator = '\n')
	writer.writerow([firstCut, secondCut])
	return outFile

def main (chrom, N, listeNucl, listeProb, distribWithLimit) :
	#Echantillonner nucléotide en accord avec leur probabilité respective
	ChoiceNuc = SampleElem (array = listeNucl, probs = listeProb)
	for j in ChoiceNuc : 
		#Echantillonner chromosome
		ChoiceChrom = random.choice (chrom)
		chromosome = ChoiceChrom.split("\t")[0]
		LenStart = ChoiceChrom.split("\t")[1]
		LenEnd = ChoiceChrom.split("\t")[2]	
		#Boucle afin de trouver la position associée au nucléotide échantillonné précédemment
		Suiteprog = False 
		while Suiteprog == False : 
			#Choix d'une positon au hasard
			ChoicePosition = random.choice (range(int(LenStart),int(LenEnd)))
			#utilisation de bedtools Getfasta pour retrouver le nucléotide associé à cette position
			Seq = createRangePybedtools(element = ChoicePosition ,chrom = chromosome)
			b = pybedtools.BedTool('''{}'''.format(Seq), from_string=True)
			fasta = pybedtools.example_filename('/groups/a2e/TAIR10/TAIR10_allChr.fasta')
			c = b.sequence(fi=fasta)
			parse = open(c.seqfn).read().split("\n")
			#print(parse)
			Getnuc = parse[1]
			for h in Getnuc : 
				#print(h)
				if h == j :
					Suiteprog = True
					pybedtools.helpers.cleanup(verbose=False,remove_all=False)
					#on passe au second cut
					#elif h == "N" : 
					#ChoicePosition = "NA"
					#Suiteprog = True
				else : 
					Suiteprog = False
					#on ré-échantillonne une position
					pybedtools.helpers.cleanup(verbose=False,remove_all=False)	
	firstCuts = ChoicePosition 
	#On obtient le second cut
	ChoiceNuc2 = SampleElem (array = listeNucl, probs = listeProb)
	for i in ChoiceNuc2 : 
		Suiteprog2 = False
		while Suiteprog2 == False : 
			#On échantillone une longueur qui suit une lois de JohnsonSB
			fragLens = fragLen (r = distribWithLimit)
			secondCuts = getsecondCut (LenEnd = int(LenEnd), firstCut = firstCuts, fragLens = int(fragLens))
			Seq2 = createRangePybedtools(element = secondCuts ,chrom = chromosome)
			b2 = pybedtools.BedTool('''{}'''.format(Seq2), from_string=True)
			fasta = pybedtools.example_filename('/groups/a2e/TAIR10/TAIR10_allChr.fasta')
			c2 = b2.sequence(fi=fasta)
			parse2 = open(c2.seqfn).read().split("\n")
			Getnuc2 = parse2[1] 
			if Getnuc2 == i :
				Suiteprog2 = True
				pybedtools.helpers.cleanup(verbose=False,remove_all=False)
			else : 
				Suiteprog2 = False
				#on ré-échantillonne une longueur
				pybedtools.helpers.cleanup(verbose=False,remove_all=False)	
	#print(secondCuts) 
	return  chromosome, firstCuts , secondCuts 

#######################################################################
################## END FUNCTIONS ######################################
#######################################################################

if __name__ == '__main__':

	##========== args
	print ('#=== Reading arguments..')
	N = int(sys.argv[1]) #Correspond au nombre de fragment voulut
	listeChromPos = Lenchr(LenchrFile = sys.argv[2])
	listeNucl, listeProb = ProbNucl(ProbNuclFile = sys.argv[3])
	outFile = sys.argv[4]
	#print(CountCut)
	#print(Countreads)
	cores=10
	#Réaliser la création de 100000 fragments suivant une loi de JohnsonSB
	distrib = JohnsonSBdistrib (a = 6 , b = 1.3)
	a = 0
	while a < 10 : 
	#Paralléliser pour faire tous les chromosomes en même temps
		result = pd.DataFrame(Parallel(n_jobs = cores)(delayed(main) (chrom = listeChromPos, N = N, listeNucl = listeNucl, listeProb=listeProb , distribWithLimit = distrib) for i in tqdm(np.arange(N))))
		sortie = open ('{}'.format(outFile), 'a')
		result.to_csv(sortie,index=False, header=False,sep="\t", doublequote= False)
		sortie.close()
		a = a + 1

