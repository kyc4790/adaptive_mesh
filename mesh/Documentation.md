# mesh/

## shapes/

Contains various ```.stl``` files to be loaded into the mdoel.

## Mesh conversion

### assignModel.m

Builds a model from an FEMesh.

### buildMesh.m

Converts a mesh from ```model.Mesh``` to the type accepted by both the constraint solver and MBO (an MBOMesh).

### getModel.m 

Builds an FEMesh from a gmsh-exported MATLAB file 

## Normal Vectors

These are solved for in getNormalVectorsFine and getNormalVectorsCoarse. Currently obtained by selecting the normal vector of an adjacent face at each control point.

## getTetDataRT

Gets some data from a tet. Called in edgeSplitXT.m
