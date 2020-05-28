# visualize/

## compareMesh.m

Takes in an old mesh and a newer, split version of the same mesh, and shows tetrahedra that have at least 2 corners that are new vertices (the result of at least two splits) by default. The number of corners can be changed by the third argument.

## overlay.m

Overlays the result of compareMesh over the result of arff's VisualizeResult, allowing the split points to be viewed simultaneously with the singular curves.

## visualizeMesh.m

Takes in a set of nodes and elements, and displays the mesh.
