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
	#{ecotype=[pathName, ...]; ecotype=[...]; ...}#
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
	#returns a set of string type nodeIDs#
	coreSet=set([])
	segmentDict=GFAfile.get_segmentDict()
	for node in segmentDict:
		if segmentDict[node].get_ecotypeNumber()>=coreNumber:
			coreSet.add(node)
	return coreSet


def get_pathTraversals(GFAfile, coreSet, coreNumber, ecotypeNumber):
	#builds a structure of all non-core bubbles and its traversals through the graph by traversing each path#
	#A bubble contains all sequence containing traversals and can have subbubbles to describe variation in the traverals#
	pathDict=GFAfile.get_pathDict()
	segmentDict=GFAfile.get_segmentDict()
	bubbleNumber=1
	for pathName in pathDict:
		pathList=GFAfile.get_path(pathName).get_pathList()
		traversal=[]
		leftAnchor=None
		leftPosition=0
		rightAnchor=None
		existingBubble=None
		pathPosition=0
		for i in range(0,len(pathList)):
			if pathList[i][:-1] in coreSet:
				if traversal:
					rightAnchor=pathList[i]
					rightPosition=pathPosition
					bubbleSet=get_traversed_bubbles(traversal, segmentDict)
					bubbleNumber, GFAfile=create_bubble(leftAnchor, leftPosition, rightAnchor, rightPosition, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile, bubbleSet)
					traversal=[]
				rightAnchor=None
				leftAnchor=pathList[i]
				leftPosition=pathPosition+GFAfile.get_segment(pathList[i][:-1]).get_sequence_length()
			else:
				traversal.append(pathList[i])
			pathPosition+=GFAfile.get_segment(pathList[i][:-1]).get_sequence_length()
		if traversal:
			bubbleSet=get_traversed_bubbles(traversal, segmentDict)
			bubbleNumber, GFAfile=create_bubble(leftAnchor, leftPosition, None, pathPosition, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile, bubbleSet)
	return GFAfile


def create_bubble(leftAnchor, leftPosition, rightAnchor, rightPosition, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile, bubbleSet):
	#adds a new traversal to the structure. If neccessary a new bubble, or subbubble is created#
	if is_novelTraversal(segmentDict, leftAnchor, rightAnchor, coreNumber):
		if leftAnchor:
			leftAnchor=segmentDict[leftAnchor[:-1]]
		if rightAnchor:
			rightAnchor=segmentDict[rightAnchor[:-1]]
		bubble=find_bubble(set(traversal), leftAnchor, rightAnchor, bubbleSet)
		if bubble:
			if bubble.get_leftAnchor()==leftAnchor and bubble.get_rightAnchor()==rightAnchor:
				bubble.add_traversal(pathName, traversal, leftPosition, rightPosition)
			else:
				subBubble=bubble.find_subBubble('X', leftAnchor, rightAnchor, set(traversal), coreNumber)
				if subBubble:
					subBubble.add_traversal(pathName, traversal, leftPosition, rightPosition)
				else:
					bubbleID=modify_bubbleID(bubble, coreNumber, ecotypeNumber)
					subBubble=GFAfile.add_bubble(bubbleID, leftAnchor, rightAnchor, set(traversal), coreNumber)
					bubble.add_subBubble(subBubble)
					subBubble.add_traversal(pathName, traversal, leftPosition, rightPosition)
					add_bubble(subBubble, traversal, segmentDict)
		else:
			bubbleNumber, bubbleID=get_bubbleID(bubbleNumber, coreNumber, ecotypeNumber)
			bubble=GFAfile.add_bubble(bubbleID, leftAnchor, rightAnchor, set(traversal), coreNumber)
			bubble.add_traversal(pathName, traversal, leftPosition, rightPosition)
			add_bubble(bubble, traversal, segmentDict)
	return bubbleNumber, GFAfile


def find_bubble(segmentSet, leftAnchor, rightAnchor, bubbleSet):
	correctBubble=False
	for bubble in bubbleSet:
		if segmentSet.issubset(bubble.get_segmentSet()):
			if correctBubble:
				if len(correctBubble.get_segmentSet())>len(bubble.get_segmentSet()):
					correctBubble=bubble
			else:
				correctBubble=bubble
		elif leftAnchor==bubble.get_leftAnchor() and rightAnchor==bubble.get_rightAnchor():
			correctBubble=bubble
			break
	return correctBubble


def add_bubble(bubble, traversal, segmentDict):
	for segment in traversal:
		segmentDict[segment[:-1]].add_bubble(bubble)
	return None


def get_traversed_bubbles(traversal, segmentDict):
	bubbleSet=set([])
	for segment in traversal:
		bubbleSet.update(segmentDict[segment[:-1]].get_bubbleList())
	return bubbleSet


def build_output(GFAfile):
	#builds the fasta-like output file#
	#>bubbleID_pathName_leftAnchor,[traversalNodes],rightAnchor#
	#concatenated fasta sequence#
	outFasta=[]
	segmentDict=GFAfile.get_segmentDict()
	bubbleList=GFAfile.get_bubbleList()
	for bubble in bubbleList:
		if bubble.get_leftAnchor():
			leftAnchor=bubble.get_leftAnchor().get_id()
		else:
			leftAnchor='None'
		if bubble.get_rightAnchor():
			rightAnchor=bubble.get_rightAnchor().get_id()
		else:
			rightAnchor='None'
		for traversal in bubble.get_traversalList():
			traversalSequence=build_traversalSequence(traversal.get_segmentList(), segmentDict)
			for path in traversal.get_pathList():
				outFasta.append('>'+bubble.get_bubbleID()+'-'+path[0]+'|'+str(path[1])+':'+str(path[2])+'|'+leftAnchor+','+','.join(traversal.get_segmentList())+','+rightAnchor)
				outFasta.append(traversalSequence)
	return outFasta


def build_traversalSequence(traversalList, segmentDict):
	traversalSequence=''
	for segment in traversalList:
		if segment[-1]=='+':
			traversalSequence+=segmentDict[segment[:-1]].get_sequence()
		elif segment[-1]=='-':
			traversalSequence+=reverseComplement(segmentDict[segment[:-1]].get_sequence())
	return traversalSequence


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


def is_novelTraversal(segmentDict, leftAnchor, rightAnchor, coreNumber):
	#checks if a traversal is a artefact with a higher core number that has already been added to the data#
	novelTraversal=False
	if leftAnchor:
		if segmentDict[leftAnchor[:-1]].get_ecotypeNumber()==coreNumber:
			novelTraversal=True
	if rightAnchor:
		if segmentDict[rightAnchor[:-1]].get_ecotypeNumber()==coreNumber:
			novelTraversal=True
	return novelTraversal


def get_bubbleID(bubbleNumber, coreNumber, ecotypeNumber):
	#builds a new bubbleID for a top-order superbubble#
	bubbleID=['0']*(ecotypeNumber-1)
	bubbleID[0]=str(bubbleNumber)
	bubbleNumber+=1
	return bubbleNumber, '.'.join(bubbleID)

		
def modify_bubbleID(bubble, coreNumber, ecotypeNumber):
	#creates a bubbleID for subbubbles by setting the correct position of the ID string to the number of current subbubbles of the parent bubble +1#
	bubbleID=bubble.get_bubbleID().split('.')
	bubbleID[ecotypeNumber-coreNumber]=str(len(bubble.get_subBubbles())+1)
	return '.'.join(bubbleID)



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
print(str(len(ecotypeDict))+' different ecotypes detected')
print('starting SV detection......')
for coreNumber in range(len(ecotypeDict), 1, -1):
	print('Current core number: '+str(coreNumber))	
	coreSet=find_coreNodes(GFAfile, ecotypeDict, coreNumber)
	GFAfile=get_pathTraversals(GFAfile, coreSet, coreNumber, len(ecotypeDict))
print('SV detection done! Constructing output file......')
outFasta=build_output(GFAfile)
write_file(outPath, '\n'.join(outFasta))
