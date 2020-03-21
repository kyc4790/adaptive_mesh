function visualizeMesh(nodes, elements)
maxVertex = max(elements(1:4, :), [], 'all');
temp = size(nodes);
maxNodes = temp(1, 2);

temp = size(elements);
numTets = temp(1,2);

edgeIndex = [1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4]';
edges = elements(edgeIndex, 1:numTets);
flattened = reshape(edges, [], 1);
edges = reshape(flattened, 2, [])';
edges = [edges edges(:, 1)];
displacement = 0;

% FEM = assembleFEMatrices(model);
% numEigenvectors = 4;
% [eigvectors, eigvalues] = eigs(FEM.K, 10, 'sm');
% colors = eigvectors(:, 10);
% colors = colors - min(colors) + eps;

a = [1:maxVertex]'; b = num2str(a); c = cellstr(b);
x = nodes(1, 1:maxVertex)';
y = nodes(2, 1:maxVertex)';
z = nodes(3, 1:maxVertex)';
tet=1:121;
scatter3(x, y, z, 100, 'filled');
hold on
% fill3(reshape(nodes(1, elements(1:4, tet)), [], 4)', reshape(nodes(2, elements(1:4, tet)), [], 4)', reshape(nodes(3, elements(1:4, tet)), [], 4)', tet);
axis equal
axis off
set(gcf, 'color', 'w');
% scatter3(x(maxVertex+1:maxNodes)', y(maxVertex+1:maxNodes)', z(maxVertex+1:maxNodes)');
text(x+displacement, y+displacement, z+ displacement, c);
patch('faces', edges, 'vertices', nodes');
hold off

end