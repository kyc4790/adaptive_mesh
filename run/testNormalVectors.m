clear; close all; clc

model = createpde('thermal', 'transient');

importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmin', 5, 'GeometricOrder', 'quadratic');
thermalProperties(model,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);

elements = model.Mesh.Elements;
nodes = model.Mesh.Nodes;
nodes = nodes - mean(nodes, 2);

maxVertex = max(elements(1:4, :), [], 'all');
temp = size(nodes);
maxNodes = temp(1, 2);
displacement = 0;

[withNormals, bdryIdx, intIdx, bdry] = getNormalVectorsCoarse(model.Mesh);

FEM = assembleFEMatrices(model);
numEigenvectors = 4;
[eigvectors, eigvalues] = eigs(FEM.K, 10, 'sm');
colors = eigvectors(:, 10);
colors = colors - min(colors) + eps;

a = [1:maxNodes]'; b = num2str(a); c = cellstr(b);
x = nodes(1, :)';
y = nodes(2, :)';
z = nodes(3, :)';
scatter3(x, y, z, 100, colors, 'filled');
hold on
axis equal
axis off
set(gcf, 'color', 'w');
% scatter3(x(maxVertex+1:maxNodes)', y(maxVertex+1:maxNodes)', z(maxVertex+1:maxNodes)');

% text(x+displacement, y+displacement, z+ displacement, c);
quiver3(nodes(1, 1:maxVertex)', nodes(2, 1:maxVertex)', nodes(3, 1:maxVertex)', withNormals(:, 1), withNormals(:, 2), withNormals(:, 3));
% patch('faces', edges, 'vertices', nodes');
hold off
