# adaptive_mesh

## External dependencies
- [arff] (https://github.com/dpa1mer/arff)

## Installation
Run ```install.m```.

## Usage

### Load a mesh. 
Load a model from a ```.stl``` file.  For example,
```matlab
    model = createpde('thermal', 'transient');
    importGeometry(model, 'sphere.stl');
    generateMesh(model, 'Hmax', 10, 'GeometricOrder', 'quadratic');
```
loads a model from ```sphere.stl``` with maximum element size 10. The ```Hmax``` argument controls how fine or coarse the mesh is.

### Adaptive Remesh
Run ```adaptiveRemesh``` on your mesh.
```matlab
    [q, model, info] = adaptiveRemesh(model)
```
For more information on the runs, run
```matlab
    [q, model, info] = adaptiveRemesh(model, [], true)
```

As an example, look at run/runOptimization.m

### Pushing Singularities Outward

#### Weighting the L and M matrices
Build a weighted mesh by running 
```matlab
    weightedMesh = buildWeightedMesh(model)
```

Adjust the weights as needed by calling
```matlab
    adjustWeights(model.Mesh, weightedMesh, @x x.^2);
```
for example, to square the weights. Then, run normal MBO.

For more information, look at run/runVariation.m

#### Relaxed Boundary MBO
Call MBO as follows:
```matlab
[q, ~, info] = MBO(mesh, OctaMBO, [], 1, 0, true, 1000, alpha);
```
where ```alpha``` is the parameter in the relaxed boundary MBO. (An ```alpha``` of zero would correspond to the unrelaxed version.)

### Visualization

#### Visualizing a Mesh from a Model

Use ```visualizeMesh``` to visualize how a single mesh looks, e.g.
```matlab
    visualizeMesh(model.Mesh.Nodes, model.Mesh.Elements)
```

#### Comparing Meshes from Two Models

Use ```compareMesh``` to highlight tetrahedra that were present in later splits but not earlier ones. E.g.
```matlab
    compareMesh(model, newModel)
```
where newModel is a newer version of model with more splits. 

Alternatively, to compare the first and last models of a split, alongside the singular curves, run
```matlab
    overlay(info)
```
Note this only works if ```adaptiveRemesh``` was called with the third argument set to ```true```.


#### Visualizing a Sequence of Adaptive Remeshes

Call ```VisualizeSequence(info)```

#### Visualize a Sequence of 

Set a variable ```all_info``` such that ```all_info(j).i``` is equal to the ```info``` output from MBO with saveIterates as on. Then, call ```visualizeVariation(all_info)```.
