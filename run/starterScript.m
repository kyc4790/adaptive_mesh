model = createpde('thermal', 'transient');
importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmin',35, 'GeometricOrder', 'quadratic');

mesh = buildMesh(model.Mesh);
qOcta = MBO(mesh, OctaMBO, [], 1, 0);