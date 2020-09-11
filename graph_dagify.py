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


def dagify_graph(GFAfile, maxSize):
	DAGset=set([])
	segmentDict=GFAfile.get_segmentDict()
	maxNodeID=get_maxNodeID(segmentDict)
	for pathName in GFAfile.get_pathDict():
		maxNodeID, DAGset, GFAfile=dagify_path(GFAfile, pathName, maxSize, maxNodeID)
	return GFAfile


def get_maxNodeID(segmentDict):
	nodeList=[]
	for nodeID in segmentDict:
		nodeList.append(int(nodeID))
	return sorted(nodeList)[-1]


def dagify_path(GFAfile, pathName, maxSize, maxNodeID):
	DAGdict={}
	DAGset=set([])
	segmentDict=GFAfile.get_segmentDict()
	pathList=GFAfile.get_path(pathName).get_pathList()
	pathPosition=0
	for i in range(len(pathList)):
		nodeID=pathList[i][:-1]
		if GFAfile.get_segmentDict()[pathList[i][:-1]].get_sequence_length()<=maxSize:
			if GFAfile.get_segmentDict()[pathList[i][:-1]].has_cycle(pathName):
				if nodeID in DAGdict:
					DAGdict[nodeID].append(pathPosition)
				else:
					DAGdict[nodeID]=[pathPosition]
					DAGset.add(nodeID)
				maxNodeID+=1
				pathList=build_new_node(GFAfile, pathName, pathPosition, pathList, i, maxNodeID)
		pathPosition+=GFAfile.get_segmentDict()[nodeID].get_sequence_length()
	clean_pathDicts(GFAfile, DAGdict, pathName)
	return maxNodeID, DAGset, GFAfile


def build_new_node(GFAfile, pathName, pathPosition, pathList, i, maxNodeID):
	newSegment=GFAfile.add_segment('\t'.join(['S', str(maxNodeID), GFAfile.get_segmentDict()[pathList[i][:-1]].get_sequence()]))
	newSegment.fill_pathDict(pathName, pathPosition)
	GFAfile.add_path(GFAfile.get_segmentDict()[pathList[i-1][:-1]], pathList[i-1][-1], newSegment, pathList[i][-1], '0M')
	if i+1<len(pathList):
		GFAfile.add_path(newSegment, pathList[i][-1], GFAfile.get_segmentDict()[pathList[i+1][:-1]], pathList[i+1][-1], '0M')
		pathList[i]=str(maxNodeID)+pathList[i][-1]
	return pathList


def clean_pathDicts(GFAfile, DAGdict, pathName):
	for segmentID in DAGdict:
		for position in DAGdict[segmentID]:
			GFAfile.get_segmentDict()[segmentID].remove_from_pathDict(pathName, position)


help="""
-h	Prints help message
-g <GFAfile>	gfa file
-r <string>	Reference sequence base (separated by '_')
-o <Filepath>	Base path and name for the generated output file
-d <int>	Dagify nodes with at max this size [default: False]
"""

DAGset=set([])
DAGsize=1
refBase=None

try:
	opts, args = getopt.getopt(sys.argv[1:],'hg:r:o:d:',['GFApath=', 'refBase=', 'DAGsize', 'outPath='])
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
	elif opt in ('-d'):
		DAGsize=int(arg)

GFAfile=gfaHandler(open_file(GFApath))
print('Graph read')
GFAfile=dagify_graph(GFAfile, DAGsize)
GFAfile=gfaHandler('\n'.join(GFAfile.rebuild_gfa()))
print('GFA dagified')
write_file(outPath, '\n'.join(GFAfile.build_gfa('\tCL:Z:changedGFA_base:'+GFApath+'_DAGsize:'+str(DAGsize)+'_refBase:'+str(refBase))))
