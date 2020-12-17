# visualize/

## Mesh Visualization

### getTetQuality.m

Gets the tet quality of all of the elements.

### visualizeMesh.m

Takes in a set of nodes and elements, and displays the mesh.

## Tools to visualize edge splitting

### compareMesh.m

Takes in an old mesh and a newer, split version of the same mesh, and shows tetrahedra that have at least 2 corners that are new vertices (the result of at least two splits) by default. The number of corners can be changed by the third argument.

### overlay.m

Overlays the result of compareMesh over the result of arff's VisualizeResult, allowing the split points to be viewed simultaneously with the singular curves.

### visualizeEdges.m

Visualizes the various costs of the edges with color. Useful in verifying whether costs are computed accurately.

## VisualizeSequence.m

Visualize a sequence of adaptive divisions. Takes in info

## VisualizeVariation.m

Plays a movie of the singular curves over time.
