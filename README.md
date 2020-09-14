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
	>1.1.0.0.0.1_AT1741_Chr1|startPosition:stopPosition|a1,nodeID+,..,..,a2

