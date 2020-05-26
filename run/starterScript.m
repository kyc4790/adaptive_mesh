model = createpde('thermal', 'transient');
importGeometry(model, 'cylinder.stl');
generateMesh(model, 'Hmax',5, 'GeometricOrder', 'quadratic');

mesh = buildMesh(model.Mesh);
qOcta = MBO(mesh, OctaMBO, [], 1, 0);