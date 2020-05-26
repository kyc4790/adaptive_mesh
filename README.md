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
