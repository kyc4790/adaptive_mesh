% model = createpde('thermal', 'transient');
% importGeometry(model, 'sphere.stl');
% generateMesh(model, 'Hmax', 10, 'GeometricOrder', 'quadratic');

[qOcta, newModel, info] = adaptiveRemesh(model, qOcta, true, 1e-3, info);
newMesh = buildMesh(newModel.Mesh);

VisualizeResult(newMesh, qOcta);
compareMesh(model, newModel);
% visualizeMesh(newModel.Mesh.Nodes, newModel.Mesh.Elements);
