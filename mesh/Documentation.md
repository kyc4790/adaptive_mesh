# mesh/

## shapes/

Contains various ```.stl``` files to be loaded into the mdoel.

## buildMesh.m

Converts a mesh from ```model.Mesh``` to the type accepted by both the constraint solver and MBO.

## Normal Vectors

These are solved for in getNormalVectorsFine and getNormalVectorsCoarse. Currently obtained by selecting the normal vector of an adjacent face at each control point.

## getTetDataRT

Gets some data from a tet. Called in edgeSplitXT.m
