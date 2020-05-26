# Documentation for adaptive_mesh

## deriv/
Contains all the files pertaining to calculating the adaptive constraints

## MBO/
Contains a copy of the MBO file from arff, without the printing.

## mesh/
Contains all of the mesh ```.stl``` files as well as tools for building the mesh.

## run/
Contains sample run files, and a few testing file.

## visualize/
Contains files for visualization.

## adaptiveRemesh.m

Performs the adaptive remesh. qOcta can either be ```[]``` for no initialization, or a coarse qOcta that corresponds to the model. saveIterates is a boolean for whether or not in-depth info is needed, and addInfo is present in case one wishes to append to an existing info struct. It returns the new info struct, the new qOcta, and the new model.

## edgeSplitXT.m

Performs the edge split, given a ```model```, its associated ```costs```, and a ```qOcta```. ```model``` is expected to be a valid model, as described in README.md. ```costs``` is expected to have length corresponding to the number of fine control points (both vertices and edges), although the points at the vertex will always have a cost of zero, such that the ```i```th entry corresponds to the control point with index ```i``` in the model. ```qOcta``` should be the coarse value of q.
