PanSV is a tool to iteratively extract non-core sequence from genome graphs.

Submodule:
- gfautils:
Install: python setup.py install

How it works:

Detecting non-core sequences by walking along each path and detecting nodes that 
are not traversed by enough ecotypes.


Class: Bubble

bubbleID: <str>
leftAnchor: <segmentObject>
rightAnchor: <segmentObject>
segmentSet: {segmentIDs} 
coreNumber: <int> (first core size for which the bubble has been detected)
parent: <bubbleObject> (None if top bubble)
subBubbles: [bubbleObjects]
