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
	coreSet=set([])
	segmentDict=GFAfile.get_segmentDict()
	for node in segmentDict:
		if segmentDict[node].get_ecotypeNumber()>=coreNumber:
			coreSet.add(node)
	return coreSet


def get_pathTraversals(GFAfile, coreSet, coreNumber, ecotypeNumber):
	#builds a list of all non-core sequences by traversing each path in a fasta format#
	pathDict=GFAfile.get_pathDict()
	segmentDict=GFAfile.get_segmentDict()
	bubbleNumber=1
	for pathName in pathDict:
		pathList=GFAfile.get_path(pathName).get_pathList()
		traversal=[]
		leftAnchor=None
		rightAnchor=None
		existingBubble=None
		for i in range(0,len(pathList)):
			if pathList[i][:-1] in coreSet:
				if traversal:
					rightAnchor=pathList[i]
					bubbleNumber, GFAfile=add_bubble(leftAnchor, rightAnchor, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile)
					traversal=[]
				rightAnchor=None
				leftAnchor=pathList[i]
			else:
				traversal.append(pathList[i])
		if traversal:
			bubbleNumber, GFAfile=add_bubble(leftAnchor, None, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile)
	for bubble in GFAfile.get_bubbleList():
		print('=======')
		print(bubble)
		print(bubble.get_bubbleID())
		if bubble.get_leftAnchor():
			print(bubble.get_leftAnchor().get_id())
		else:
			print('None')
		if bubble.get_rightAnchor():
			print(bubble.get_rightAnchor().get_id())
		else:
			print('None')
		print(bubble.get_subBubbles())
	return GFAfile


def add_bubble(leftAnchor, rightAnchor, traversal, pathName, segmentDict, coreNumber, ecotypeNumber, bubbleNumber, GFAfile):
	if is_novelTraversal(segmentDict, leftAnchor, rightAnchor, coreNumber):
		if leftAnchor:
			leftAnchor=segmentDict[leftAnchor[:-1]]
		if rightAnchor:
			rightAnchor=segmentDict[rightAnchor[:-1]]
		bubble=GFAfile.has_bubble(traversal, leftAnchor, rightAnchor)
		if bubble:
			if bubble.get_leftAnchor()==leftAnchor and bubble.get_rightAnchor()==rightAnchor:
				bubble.add_traversal(pathName, traversal)
			else:
				subBubble=bubble.find_subBubble('X', leftAnchor, rightAnchor, set(traversal), coreNumber)
				if subBubble:
					subBubble.add_traversal(pathName, traversal)
				else:
					bubbleID=modify_bubbleID(bubble, coreNumber, ecotypeNumber)
					subBubble=GFAfile.add_bubble(bubbleID, leftAnchor, rightAnchor, set(traversal), coreNumber)
					bubble.add_subBubble(subBubble)
					subBubble.add_traversal(pathName, traversal)
		else:
			bubbleNumber, bubbleID=get_bubbleID(bubbleNumber, coreNumber, ecotypeNumber)
			bubble=GFAfile.add_bubble(bubbleID, leftAnchor, rightAnchor, set(traversal), coreNumber)
			bubble.add_traversal(pathName, traversal)
	return bubbleNumber, GFAfile

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
	novelTraversal=False
	if leftAnchor:
		if segmentDict[leftAnchor[:-1]].get_ecotypeNumber()==coreNumber:
			novelTraversal=True
	if rightAnchor:
		if segmentDict[rightAnchor[:-1]].get_ecotypeNumber()==coreNumber:
			novelTraversal=True
	return novelTraversal


def get_bubbleID(bubbleNumber, coreNumber, ecotypeNumber):
	bubbleID=['0']*(ecotypeNumber-1)
	bubbleID[0]=str(bubbleNumber)
	bubbleNumber+=1
	return bubbleNumber, '.'.join(bubbleID)

		
def modify_bubbleID(bubble, coreNumber, ecotypeNumber):
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
for coreNumber in range(len(ecotypeDict), 1, -1):
	coreSet=find_coreNodes(GFAfile, ecotypeDict, coreNumber)
	GFAfile=get_pathTraversals(GFAfile, coreSet, coreNumber, len(ecotypeDict))
#outFasta=build_output(GFAfile)
#write_file(outPath, '\n'.join(traversalList))
