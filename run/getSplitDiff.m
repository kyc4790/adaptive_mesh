model = createpde('thermal', 'transient');
importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmax', 3, 'GeometricOrder', 'quadratic');

mesh = buildMesh(model.Mesh);
[qOcta, ~, info] = MBO(mesh, OctaMBO);

IHh = buildIHh(model.Mesh);
qFine = (IHh*qOcta')';

[newModel, newQ] = edgeSplitXT(model, randperm(max(model.Mesh.Elements,[], 'all'))', qFine, mesh, 0.05);
newMesh = buildMesh(newModel.Mesh); 
[qOctaNew, ~, infoNew] = MBO(newMesh, OctaMBO, newQ);