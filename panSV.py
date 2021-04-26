#!/usr/bin/env python

import sys
import getopt
import argparse
from gfa import gfaHandler


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


def get_pathTraversals(GFAfile, ecotypeNumber):
	#builds a structure of all non-core bubbles and its traversals through the graph by traversing each path#
	#A bubble contains all sequence containing traversals and can have subbubbles to describe variation in the traverals#
	pathDict=GFAfile.get_pathDict()
	segmentDict=GFAfile.get_segmentDict()
	for pathName in pathDict:
		openTraversalDict=initializeTraversalDict(ecotypeNumber)
		closedTraversalDict=initializeTraversalDict(ecotypeNumber)
		pathList=GFAfile.get_path(pathName).get_pathList()
		leftAnchor=None
		rightAnchor=None
		pathPosition=0
		if segmentDict[pathList[0][:-1]].get_ecotypeNumber()!=ecotypeNumber:
			openTraversalDict[ecotypeNumber]=[None, None, 0, None, []]
		for i in range(0,len(pathList)-1):
			pathPosition+=segmentDict[pathList[i][:-1]].get_sequence_length()
			if segmentDict[pathList[i][:-1]].get_ecotypeNumber()!=segmentDict[pathList[i+1][:-1]].get_ecotypeNumber():
				if segmentDict[pathList[i][:-1]].get_ecotypeNumber()<segmentDict[pathList[i+1][:-1]].get_ecotypeNumber():
					openTraversalDict=add_segment(openTraversalDict, pathList[i])
					# close all bubbles of lower degree #
					openTraversalDict, closedTraversalDict, GFAfile=closeTraversal(openTraversalDict, segmentDict[pathList[i+1][:-1]], pathPosition, closedTraversalDict, ecotypeNumber, GFAfile, pathName)
				elif segmentDict[pathList[i][:-1]].get_ecotypeNumber()>segmentDict[pathList[i+1][:-1]].get_ecotypeNumber():
					openTraversalDict=add_segment(openTraversalDict, pathList[i])
					# open new bubble #
					openTraversalDict[segmentDict[pathList[i][:-1]].get_ecotypeNumber()]=[segmentDict[pathList[i][:-1]], None, pathPosition, None, []]
			else:
				if i!=0:
					# add node to all open bubbles #
					openTraversalDict=add_segment(openTraversalDict, pathList[i])
		if segmentDict[pathList[-1][:-1]].get_ecotypeNumber()!=ecotypeNumber:
			openTraversalDict=add_segment(openTraversalDict, pathList[-1])
		openTraversalDict, closedTraversalDict, GFAfile=closeTraversal(openTraversalDict, segmentDict[pathList[i+1][:-1]], pathPosition, closedTraversalDict, ecotypeNumber, GFAfile, pathName)
	# add unique IDs to bubbles #
	GFAfile=nameBubbles(GFAfile, ecotypeNumber)
	return GFAfile


def initializeTraversalDict(ecotypeNumber):
	#{coreLevel=[leftAnchor, rightAnchor, startPosition, stopPosition, segmentList],..}
	traversalDict={}
	for i in range(2, ecotypeNumber):
		traversalDict[i]=[]
	traversalDict[ecotypeNumber]=[]
	return traversalDict


def add_segment(traversalDict, segment):
	for coreLevel in traversalDict:
		if traversalDict[coreLevel]:
			traversalDict[coreLevel][4].append(segment)
	return traversalDict


def closeTraversal(traversalDict, segmentObject, pathPosition, closedTraversalDict, ecotypeNumber, GFAfile, pathName):
	print(traversalDict)
	print(segmentObject)
	for coreLevel in traversalDict:
		if coreLevel<=segmentObject.get_ecotypeNumber():
			if traversalDict[coreLevel]:
				traversalDict[coreLevel][1]=segmentObject
				traversalDict[coreLevel][3]=pathPosition
				closedTraversalDict[coreLevel].append(traversalDict[coreLevel])
				if segmentObject.get_ecotypeNumber()==ecotypeNumber and coreLevel==ecotypeNumber:
					closedTraversalDict, GFAfile=create_bubbles(closedTraversalDict, GFAfile, pathName, ecotypeNumber)
				traversalDict[coreLevel]=[]
			else:
				if coreLevel==segmentObject.get_ecotypeNumber() and coreLevel!=ecotypeNumber:
					traversal=False
					traversalSearchLevel=coreLevel
					while not traversal:
						traversalSearchLevel+=1
						if ecotypeNumber>=traversalSearchLevel:
							if traversalDict[traversalSearchLevel]:
								traversal=traversalDict[traversalSearchLevel].copy()
								traversal[1]=segmentObject
								traversal[3]=pathPosition
								traversal[4]=traversalDict[traversalSearchLevel][4].copy()
								closedTraversalDict[segmentObject.get_ecotypeNumber()].append(traversal)
						else:
							break
	return traversalDict, closedTraversalDict, GFAfile


def create_bubbles(closedTraversalDict, GFAfile, pathName, ecotypeNumber):
	parentList=[]
	for i in range(ecotypeNumber, 1, -1):
		for traversal in closedTraversalDict[i]:
			bubble=bubble_exists(traversal[0], traversal[1])
			if bubble:
				bubble.add_traversal(pathName, traversal[4], traversal[2], traversal[3])
			else:
				parent=find_parent(parentList, traversal[4])
				# create new bubble and add traversal#
				bubble=GFAfile.add_bubble(None, traversal[0], traversal[1], traversal[4], i, parent)
				bubble.add_traversal(pathName, traversal[4], traversal[2], traversal[3])
			parentList.append(bubble)
		closedTraversalDict[i]=[]
	return closedTraversalDict, GFAfile


def bubble_exists(leftAnchor, rightAnchor):
	bubble=None
	if leftAnchor and rightAnchor:
		if set(leftAnchor.get_leftAnchor()).intersection(rightAnchor.get_rightAnchor()):
			bubble=list(set(leftAnchor.get_leftAnchor()).intersection(rightAnchor.get_rightAnchor()))[0]
	return bubble


def find_parent(parentList, traversalList):
	parent=None
	for bubble in parentList:
		if set(traversalList).issubset(bubble.get_segmentSet()):
			if bubble.get_subBubbles():
				parent=find_parent(bubble.get_subBubbles(), traversalList)
				if not parent:
					parent=bubble
				break
			else:
				parent=bubble
				break
	return parent


def nameBubbles(GFAfile, ecotypeNumber):
	bubbleDict=GFAfile.get_bubbleDict()
	nameList=[0]*(ecotypeNumber-1)
	for bubble in bubbleDict[str(ecotypeNumber)]:
		nameList=bubble.set_bubbleID(nameList)
	return GFAfile


### add functionality!!! ###
def add_bubble(bubble, traversal, segmentDict):
	for segment in traversal:
		segmentDict[segment[:-1]].add_bubble(bubble)
	return None


def get_traversed_bubbles(traversal, segmentDict):
	bubbleSet=set([])
	for segment in traversal:
		bubbleSet.update(segmentDict[segment[:-1]].get_bubbleList())
	return bubbleSet


def getPAVtraversals(GFAfile):
	pathDict=GFAfile.get_pathDict()
	for pathName in pathDict:
		pathList=GFAfile.get_path(pathName).get_pathList()
		pathPosition=0
		for i in range(0,len(pathList)):
			if GFAfile.get_segment(pathList[i][:-1]).get_rightAnchor():
#				for bubble in GFAfile.get_segment(pathList[i-1][:-1]).get_leftAnchor():
				for bubble in GFAfile.get_segment(pathList[i][:-1]).get_rightAnchor():
					if bubble in GFAfile.get_segment(pathList[i-1][:-1]).get_leftAnchor():
						leftPosition=pathPosition
						rightPosition=pathPosition+1
#						leftPosition=pathPosition+GFAfile.get_segment(pathList[i][:-1]).get_sequence_length()+1
#						rightPosition=pathPosition+GFAfile.get_segment(pathList[i][:-1]).get_sequence_length()+1
						bubble.add_traversal(pathName, 'PAV', leftPosition, rightPosition)
						break
			pathPosition+=GFAfile.get_segment(pathList[i][:-1]).get_sequence_length()
	return GFAfile


def build_output(GFAfile):
	#builds the fasta-like output file#
	#>bubbleID_pathName_leftAnchor,[traversalNodes],rightAnchor#
	#concatenated fasta sequence#
	outFasta=[]
	outBED=[]
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
			if isinstance(traversal.get_segmentList(), list) or isinstance(traversal.get_segmentList(), set):
				traversalSequence=build_traversalSequence(traversal.get_segmentList(), segmentDict)
			else:
				traversalSequence=traversal.get_segmentList()
			for path in traversal.get_pathList():
				outBED.append('\t'.join([path[0], str(path[1]), str(path[2]), bubble.get_bubbleID(), str(bubble.get_coreNumber())]))
				if isinstance(traversal.get_segmentList(), list) or isinstance(traversal.get_segmentList(), set):
					outFasta.append('>'+bubble.get_bubbleID()+'-'+path[0]+'|'+str(path[1])+':'+str(path[2])+'|'+leftAnchor+','+','.join(traversal.get_segmentList())+','+rightAnchor)
				else:
					outFasta.append('>'+bubble.get_bubbleID()+'-'+path[0]+'|'+str(path[1])+':'+str(path[2])+'|'+leftAnchor+','+rightAnchor)
				outFasta.append(traversalSequence)
	return outFasta, outBED


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


#def is_novelTraversal(segmentDict, leftAnchor, rightAnchor, coreNumber):
#	#checks if a traversal belongs to an existing bubble with a higher core number that has already been added to the data#
#	novelTraversal=False
#	if leftAnchor:
#		if segmentDict[leftAnchor[:-1]].get_ecotypeNumber()==coreNumber:
#			novelTraversal=True
#	if rightAnchor:
#		if segmentDict[rightAnchor[:-1]].get_ecotypeNumber()==coreNumber:
#			novelTraversal=True
#	return novelTraversal


def get_bubbleID(bubbleNumber, coreNumber, ecotypeNumber):
	#builds a new bubbleID for a top-order superbubble#
	bubbleID=['0']*(ecotypeNumber-1)
	bubbleID[ecotypeNumber-coreNumber]=str(bubbleNumber)
	bubbleNumber+=1
	return bubbleNumber, '.'.join(bubbleID)

		
def modify_bubbleID(bubble, coreNumber, ecotypeNumber):
	#creates a bubbleID for subbubbles by setting the correct position of the ID string to the number of current subbubbles of the parent bubble +1#
	bubbleID=bubble.get_bubbleID().split('.')
	bubbleID[ecotypeNumber-coreNumber]=str(len(bubble.get_subBubbles())+1)
	return '.'.join(bubbleID)


def get_outStats(GFAfile):
	outStats=['\t'.join(['bubbleID', 'coreNumber', 'subBubbles', 'sequence', 'minLen', 'maxLen', 'avgLen', 'traversals', 'pathTraversals'])]
	bubbleList=GFAfile.get_bubbleList()
	for bubble in bubbleList:
		bubbleID=bubble.get_bubbleID()
		coreNumber=bubble.get_coreNumber()
		subBubbles=str(len(bubble.get_subBubbles()))
		sequence=str(get_bubbleSequence(bubble, GFAfile))
		minLen, maxLen, avgLen=get_traversalLengths(bubble.get_traversalList(), GFAfile)
		traversals=str(len(bubble.get_traversalList()))
		pathTraversals=str(get_pathTraversalNumber(bubble.get_traversalList()))
		outStats.append('\t'.join([bubbleID, coreNumber, subBubbles, sequence, str(minLen), str(maxLen), str(avgLen), traversals, pathTraversals]))
	return outStats


def get_pathTraversalNumber(traversalList):
	pathTraversals=0
	for traversal in traversalList:
		pathTraversals+=len(traversal.get_pathList())
	return pathTraversals


def get_bubbleSequence(bubble, GFAfile):
	bubbleSequence=0
	for segment in bubble.get_segmentSet():
		bubbleSequence+=GFAfile.get_segment(segment[:-1]).get_sequence_length()
	return bubbleSequence


def get_traversalLengths(traversalList, GFAfile):
	minLen=None
	maxLen=0
	avgLen=0
	combLen=0
	for traversal in traversalList:
		traversalLength=0
		if isinstance(traversal.get_segmentList(), list) or isinstance(traversal.get_segmentList(), set):
			for segment in traversal.get_segmentList():
				segmentLength=GFAfile.get_segment(segment[:-1]).get_sequence_length()
				traversalLength+=segmentLength
		else:
			traversalLength=0
		combLen+=traversalLength
		if minLen:
			if traversalLength<minLen:
				minLen=traversalLength
		else:
			minLen=traversalLength
		if traversalLength>maxLen:
			maxLen=traversalLength
	avgLen=combLen/float(len(traversalList))		
	return minLen, maxLen, avgLen

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "PanSV: Pan-genome SV detection algorithm by C. Kubica")
	parser.add_argument("-g", "--gfa", help="gfa file", required=True)
	parser.add_argument("-o", "--output", help="Base path and name for the generated output file ", required=True)
	args = parser.parse_args();

	GFApath = args.gfa
	outPath = args.output



	GFAfile=gfaHandler(open_file(GFApath))
	print('Graph read')
	ecotypeDict=get_ecotypeDict(GFAfile.get_pathDict())
	print(str(len(ecotypeDict))+' different ecotypes detected')
	print('starting SV detection......')
	GFAfile=get_pathTraversals(GFAfile, len(ecotypeDict))
	print('Detecting PAV traversals...')
	GFAfile=getPAVtraversals(GFAfile)
	print('SV detection done! Constructing output files......')
	outFasta, outBED=build_output(GFAfile)
	outStats=get_outStats(GFAfile)
	write_file(outPath+'.fasta', '\n'.join(outFasta))
	write_file(outPath+'.bed', '\n'.join(outBED))
	write_file(outPath+'.stats', '\n'.join(outStats))
