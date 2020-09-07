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
	ecotypeDict={}
	for path in pathList:
		if path.split('_')[0] in ecotypeDict:
			ecotypeDict[path.split('_')[0]].append(path)
		else:
			ecotypeDict[path.split('_')[0]]=[path]
	return ecotypeDict


def find_coreNodes(GFAfile, ecotypeDict, coreNumber):
	coreSet=set([])
	segmentDict=GFAfile.get_segmentDict()
	for node in segmentDict:
		if segmentDict[node].get_ecotypeNumber()>=coreNumber:
			coreSet.add(node)
	return coreSet


def get_pathTraversals(GFAfile, coreSet):
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
					traversalList.append(pathName+'\t'+str(len(SVsequence)))
#					traversalList.append('>'+pathName+':'+';'.join(SVheader))
#					traversalList.append(SVsequence)
					SVheader=[]
					SVsequence=''
			else:
				SVheader.append(pathList[i])
				if pathList[i][-1]=='+':
					SVsequence+=segmentDict[pathList[i][:-1]].get_sequence()
				else:
					SVsequence+=reverseComplement(segmentDict[pathList[i][:-1]].get_sequence())
		if SVheader:
			traversalList.append(pathName+'\t'+str(len(SVsequence)))
#			traversalList.append('>'+pathName+':'+';'.join(SVheader))
#			traversalList.append(SVsequence)
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
			
		


def prune_path(GFAfile, pathName, maxSize, refBase):
	segmentDict=GFAfile.get_segmentDict()
	pathList=GFAfile.get_path(pathName).get_pathList()
	pathPosition=segmentDict[pathList[0][:-1]].get_sequence_length()
	changedPathPosition=pathPosition
	refLength=0
	noRefLength=0
	for i in range(0,len(pathList)-1):
		nextNode=None
		oldPath=None
		if pathList[i]!='indel':
			if pathList[i-1]!='indel' and len(pathList[i-1].split(','))==1:
				previousNode=pathList[i-1][:-1]
			nodeID=pathList[i][:-1]
			followingNode=pathList[i+1][:-1]
			if segmentDict[nodeID].has_path(refBase + '_' + pathName.split('_')[-1]):
				if segmentDict[followingNode].has_path(refBase + '_' + pathName.split('_')[-1]):
					nextNode=prune_refIndel(segmentDict[nodeID], segmentDict[followingNode], refBase+'_'+pathName.split('_')[-1], pathList[i], maxSize)
				else:
					nextNode, oldPath=prune_refSV(segmentDict[nodeID], refBase+'_'+pathName.split('_')[-1], pathList[i], pathList, i+1, segmentDict, maxSize)
			elif not is_indel(segmentDict[nodeID], segmentDict[followingNode], pathName) and is_bubble(segmentDict[previousNode], segmentDict[nodeID], segmentDict[followingNode]) and segmentDict[nodeID].get_sequence_length()<=maxSize:
				nextNode=prune_nonReference_bubble(segmentDict[previousNode], segmentDict[nodeID], segmentDict[followingNode], pathName, refBase, maxSize)
		if nextNode:
			pathList[i]=nextNode
			if oldPath:
				j=i+1
				for segment in oldPath:
					if segment.get_id()==pathList[j][:-1]:
						pathList[j]='indel'
					j+=1
	newPath, newCigar=clean_path(pathList, pathName)
	GFAfile.get_pathDict()[pathName].change_pathList(newPath, newCigar)
	return GFAfile


def prune_refIndel(currentNode, followingNode, refPath, pathString, maxSize):
	nextNode=None
	traversalList=reconstruct_traversals(currentNode, followingNode, refPath, maxSize)
	shortestTraversal=None
	traversalLength=maxSize
	for traversal in traversalList:
		if traversal[0]<traversalLength:
			shortestTraversal=traversal
			traversalLength=traversal[0]
	if shortestTraversal:
		nextNode=pathString+','+build_traversalString(currentNode, shortestTraversal[1], followingNode)
	return nextNode


def prune_refSV(leftAnchor, refPath, pathString, pathList, i, segmentDict, maxSize):
	nextNode=None
	oldTraversal=[]
	traversalLength=0
	rightAnchor=None
	closedBubble=True
	while not segmentDict[pathList[i][:-1]].has_path(refPath):
		traversalLength+=segmentDict[pathList[i][:-1]].get_sequence_length()
		oldTraversal.append(segmentDict[pathList[i][:-1]])
		i+=1
		if i>=len(pathList):
			closedBubble=False
			break
		if traversalLength>maxSize:
			closedBubble=False
			break
	if closedBubble:
		rightAnchor=segmentDict[pathList[i][:-1]]
		if traversalLength<=maxSize and i!=len(pathList):
			traversalList=reconstruct_traversals(leftAnchor, rightAnchor, refPath, maxSize)
			if traversalList:
				shortestTraversal=None
				traversalLength=maxSize
				for traversal in traversalList:
					if traversal[0]<traversalLength:
						shortestTraversal=traversal
						traversalLength=traversal[0]
				if shortestTraversal:
					nextNode=pathString+','+build_traversalString(leftAnchor, shortestTraversal[1], rightAnchor)
	return nextNode, oldTraversal


def reconstruct_traversals(leftAnchor, rightAnchor, pathName, maxSize=10000000000):
	###[[traversalLength, nodeList],[...],...]###
	traversalList=[]
	for pathPosition in leftAnchor.get_path_positions(pathName):
		traversalLength=0
		nodeList=[]
		nextNode=leftAnchor
		while nextNode!=rightAnchor:
			nextNode, nextNodePosition=nextNode.get_successorNode(pathName, pathPosition)
			if nextNode:
				traversalLength+=nextNode.get_sequence_length()
				nodeList.append(nextNode)
			if traversalLength>maxSize:
				nodeList=[]
				break
			else:
				break
		if nodeList:
			traversalList.append([traversalLength, nodeList])
	return traversalList


def build_traversalString(leftAnchor, nodeList, rightAnchor):
	traversalList=[]
	for successor in leftAnchor.get_successors():
		if successor.get_rightSegment()==nodeList[0]:
			traversalList.append(nodeList[0].get_id()+successor.get_rightOrientation())
	for i in range(len(nodeList)-1):
		for successor in nodeList[i].get_successors():
			if successor.get_rightSegment()==nodeList[i+1]:
				traversalList.append(nodeList[i+1].get_id()+successor.get_rightOrientation())
	return ','.join(traversalList)

			
def prune_nonReference_bubble(leftAnchor, bubbleNode, rightAnchor, pathName, refBase, maxSize):
	nextNode=None
	traversalList=[]
	inDel=False
	for successor in leftAnchor.get_successors():
		if successor.get_rightSegment()!=rightAnchor and successor.get_rightSegment()!=bubbleNode:
			for sndSuccessor in successor.get_rightSegment().get_successors():
				if sndSuccessor.get_rightSegment()==rightAnchor:
					if sndSuccessor.get_rightSegment().get_sequence_length()<maxSize:
						traversalList.append(successor.get_rightSegment())
		elif successor.get_rightSegment()==rightAnchor:
			inDel=True
			break
	if inDel:
		nextNode='indel'
	elif traversalList:
		traversalNumber=len(bubbleNode.get_pathDict())
		successorObject=None
		for traversal in traversalList:
			if traversalNumber<=len(traversal.get_pathDict()) and traversal.get_sequence_length()<maxSize:
				traversalNumber=len(traversal.get_pathDict())
				successorObject=traversal
		if successorObject:
			for successor in leftAnchor.get_successors():
				if successor.get_rightSegment()==successorObject:
					nextNode=successorObject.get_id()+successor.get_rightOrientation()
	return nextNode


def clean_path(pathList, pathName):
	cleanPath=[]
	for id in pathList:
		if id!='indel':
			cleanPath.extend(id.split(','))
	return cleanPath, ['0M']*len(cleanPath)


def is_indel(leftAnchor, rightAnchor, pathName):
	indel=False
	for path in leftAnchor.get_pathDict():
		if path!=pathName and rightAnchor.has_path(path):
			for pathPosition in leftAnchor.get_path_positions(path):
				if pathPosition+leftAnchor.get_sequence_length() in rightAnchor.get_path_positions(path):
					indel=True
					break
	return indel


def is_bubble(leftAnchor, bubbleNode, rightAnchor):
	bubble=False
	for successor in leftAnchor.get_successors():
		if successor.get_rightSegment()!=rightAnchor and successor.get_rightSegment()!=bubbleNode:
			for sndSuccessor in successor.get_rightSegment().get_successors():
					if sndSuccessor.get_rightSegment()==rightAnchor:
						bubble=True
						break
	return bubble







help=""" 
-h		Prints help message 
-g <GFAfile>	gfa file 
-r <string>	Reference sequence base (separated by '_') 
-c <int>	Number of ecotypes that need to be present for a node to be considered as core Prune nodes with this max size [default=2]
-o <Filepath>	Base path and name for the generated output file 
"""

corenNumber=2

try:
	opts, args = getopt.getopt(sys.argv[1:],'hg:r:c:o:',['GFApath=', 'refBase=', 'coreNumber', 'outPath='])
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
	elif opt in ('-r'):
		refBase=arg
	elif opt in ('-c'):
		coreNumber=int(arg)


GFAfile=gfaHandler(open_file(GFApath))
print('Graph read')
ecotypeDict=get_ecotypeDict(GFAfile.get_pathDict())
coreSet=find_coreNodes(GFAfile, ecotypeDict, coreNumber)
print(len(coreSet))
traversalList=get_pathTraversals(GFAfile, coreSet)
write_file(outPath, '\n'.join(traversalList))
