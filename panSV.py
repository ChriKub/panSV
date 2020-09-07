#!/usr/bin/env python


import sys
import getopt
from gfa import gfaHandler

from datetime import datetime


def open_file(inpath):
	file = open(inpath, 'r')
	input = file.read()
	file.close()
	return input


def write_file(outpath, output):
	file = open(outpath, 'w')
	file.write(output)
	file.close()
	return None


def get_ecotypeDict(pathList):
	#creates a dictionary of ecotypes by splitting the path name at _ #
	#ecotype_pathID#
	ecotypeDict={}
	for path in pathList:
		if path.split('_')[0] in ecotypeDict:
			ecotypeDict[path.split('_')[0]].append(path)
		else:
			ecotypeDict[path.split('_')[0]]=[path]
	return ecotypeDict


def find_coreNodes(GFAfile, ecotypeDict, coreNumber): 
	#builds a set that contains all nodes meeting the current core requirements#
	coreSet={}
	segmentDict=GFAfile.get_segmentDict()
	for node in segmentDict:
		if segmentDict[node].get_ecotypeNumber()>=coreNumber:
			coreSet.add(node)
	return coreSet


def get_pathTraversals(GFAfile, coreSet):
	#builds a list of all non-core sequences by traversing each path in a fasta format#
	traversalList=[]
	pathDict=GFAfile.get_pathDict()
	segmentDict=GFAfile.get_segmentDict()
	for pathName in pathDict:
		pathList=GFAfile.get_path(pathName).get_pathList()
		SVheader=[]
		SVsequence=''
		for i in range(0,len(pathList)):
			if pathList[i][:-1] in coreSet:
				if SVheader:
					traversalList.append('>'+pathName+':'+';'.join(SVheader))
					traversalList.append(SVsequence)
					SVheader=[]
					SVsequence=''
			else:
				SVheader.append(pathList[i])
				if pathList[i][-1]=='+':
					SVsequence+=segmentDict[pathList[i][:-1]].get_sequence()
				else:
					SVsequence+=reverseComplement(segmentDict[pathList[i][:-1]].get_sequence())
		if SVheader:
			traversalList.append('>'+pathName+':'+';'.join(SVheader))
			traversalList.append(SVsequence)
	return traversalList



def reverseComplement(sequence):
	reverse=''
	for base in sequence[::-1]:
		if base.upper()=='A':
			reverse+='T'
		elif base.upper()=='T':
			reverse+='A'
		elif base.upper()=='C':
			reverse+='G'
		elif base.upper()=='G':
			reverse+='C'
		else:
			reverse+='N'
	return reverse
			
		



help=""" 
-h		Prints help message 
-g <GFAfile>	gfa file 
-o <Filepath>	Base path and name for the generated output file 
"""

try:
	opts, args = getopt.getopt(sys.argv[1:],'hg:r:c:o:',['GFApath=', 'outPath='])
except getopt.GetoptError:
	print(help)
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print(help)
		sys.exit()
	elif opt in ('-g'):
		GFApath=arg
	elif opt in ('-o'):
		outPath=arg


GFAfile=gfaHandler(open_file(GFApath))
print('Graph read')
ecotypeDict=get_ecotypeDict(GFAfile.get_pathDict())
coreNumber=len(ecotypeDict)
coreSet=find_coreNodes(GFAfile, ecotypeDict, coreNumber)
traversalList=get_pathTraversals(GFAfile, coreSet)
write_file(outPath, '\n'.join(traversalList))
