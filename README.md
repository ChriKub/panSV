PanSV is a tool to iteratively extract non-core sequence from genome graphs.

Submodule:
- gfautils:
Install: python setup.py install

How it works:

Detecting non-core sequences by walking along each path and detecting nodes that 
are not traversed by enough ecotypes.


# **Class: Bubble**

* bubbleID: <str>
* leftAnchor: <segmentObject>
* rightAnchor: <segmentObject>
* segmentSet: {segmentIDs} 
* coreNumber: <int> (first core size for which the bubble has been detected)
* parent: <bubbleObject> (None if top bubble)
* subBubbles: [bubbleObjects]

## Methods:
> *get_bubbleID()*
> returns the ID string of the bubble

> *get_Anchors()*
> returns both anchors as segmentObjects: leftAnchor, rightAnchor

> *get_leftAnchor()*
> returns the leftAnchor as segmentObject

> *get_rightAnchor()*
> returns the rightAnchor as segmentObject

> *get_segmentSet()*
> returns a set of all segmentIDs that are part of this bubble

> *find_subBubble(bubbleID: <str>, leftAnchor: <segmentObject>, rightAnchor: <segmentObject>, traversalSet: <set>, coreNumber: <int>)
> Returns a existing sub bubble that fits the current traversal, if non can be found **None**

> *get_subBubbles()*
> Returns a list of all bubbleObjects that are assigned as sub bubbles.

> *add_segments(segmentList: <list>)*
> Adds all segmentIDs in the segmentList to the bubbles segmentSet

> *add_traversal(pathName: <string>, segmentList: <list>)*
> Adds a new traversalObject to this bubble, or if this traversal already exists adds the path to the existing traversal

> *get_traversalList()*
> Returns a list of all traversalObjects attached to the bubbleObject

> *add_subBubble(subBubble: <bubbleObject>)*
> Adds the new subBubble to the bubbles subBubbleList. 
