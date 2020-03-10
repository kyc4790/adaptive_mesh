model = createpde('thermal', 'transient');
importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmin', 1, 'GeometricOrder', 'quadratic');

nodes = model.Mesh.Nodes;
elements = model.Mesh.Elements;
nodes = nodes - mean(nodes, 2);

temp = size(nodes);
maxNodes = temp(1, 2);
maxVertex = max(elements(1:4, :), [], 'all');

mesh = buildMesh(model.Mesh);


qOcta = MBO(mesh, OctaMBO, [], 1, 0);
qOcta = OctaManopt(mesh, qOcta);