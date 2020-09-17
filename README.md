PanSV is a tool to iteratively extract non-core sequence from genome graphs. It 
builds a tree structure detailing the top-down variants in super-bubbles.

#How to use:

	./panSV.py -g INGFA -o OUTFASTA



##Submodule:
- gfautils:
Install: python setup.py install


#How it works:

* Detect core nodes (based on ecotype traversals)
* Iterate over all paths
  * traverse the path from left to right
    * if node is core:
      * store node order until next core node
      * search for existing bubbles with the same anchors
        * add as traversal to existing bubble
      * if no existing bubble: search for larger bubble with full node subset
        * add new sub-bubble to parent bubble
* repeat with core size -1
* write all path supported bubble traversals in a fasta-like format
  * naming: 1.0.0.0.1 
  * increment bubble, or subbubble number based on core size step and number of existing subbubles in this step
    * Fasta header (identified by anchor)

##Output Files:

PanSV produces 3 different output files. A fasta file that contains all non-core 
sequence for each path. A bed file specifying the positions of bubbles for easy 
access and a stats file giving stats for each bubble type.

###Fasta output:

The fasta header specifies the exact bubble ID and additional information on 
path, position and traversed nodes.
	>bubbleID_pathID|startPosition:stopPosition|leftAnchor,nodeList...,rightAnchor
	>1.1.0.0.4.0_AT1741_Chr1|143581:143645|1580,1581+,1582+,1583+,1584+,1585+,1586+,1587

###BED output:

The BED file has the standard bed format. The name field contains the bubbleID 
and the score field the core set size of the bubble.

###stats output:

The stats output gives basic information for each detected bubble in a tab 
separated format.

1. _bubbleID_       Unique bubble identifyer
2. _coreNumber_     Core set size
3. _subBubbles_     Number of direct sub bubbles
4. _sequence_       Concatenated length of all nodes that form this bubble
5. _minLen_         Length of the shortest traversal
6. _maxLen_         Length of the logest traversal
7. _avgLen_         Average traversal length (combined traversal length / number of traversals)
8. _traversals_     Number of unique sequence containing traversals
9. _pathTraversals_ Number of times that path traverse through sequence containing traverals
