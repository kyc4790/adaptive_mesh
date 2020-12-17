function visualizeEdges(nodes, elements, colors)
maxVertex = max(elements(1:4, :), [], 'all');
temp = size(nodes);
maxNodes = temp(1, 2);

temp = size(elements);
numTets = temp(1,2);

edgeIndexes = [1, 2, 5; 2, 3, 6; 1, 3, 7; 1, 4, 8; 2, 4, 9; 3, 4, 10]';
edges = reshape(elements(edgeIndexes, :), 3, []);
edges = unique(sort(edges)', 'rows');

% FEM = assembleFEMatrices(model);
% numEigenvectors = 4;
% [eigvectors, eigvalues] = eigs(FEM.K, 10, 'sm');
% colors = eigvectors(:, 10);
% colors = colors - min(colors) + eps;
% maxNodes = 8;
a = [1:maxNodes]'; b = num2str(a); c = cellstr(b);
x = nodes(1, maxVertex+1:maxNodes)';
x = x - mean(x);
y = nodes(2, maxVertex+1:maxNodes)';
y = y - mean(y);
z = nodes(3, maxVertex+1:maxNodes)';
z = z - mean(z);


if nargin < 3
    colors = zeros(maxNodes, 1);
end

edge_colors = colors(maxVertex+1:maxNodes);
vals= edge_colors > 0;

scatter3(x(vals), y(vals), z(vals), edge_colors(vals) / max(edge_colors(vals)) * 100, edge_colors(vals) / max(edge_colors(vals)) * 2.5e4, 'filled');

hold on
% fill3(reshape(nodes(1, elements(1:4, tet)), [], 4)', reshape(nodes(2, elements(1:4, tet)), [], 4)', reshape(nodes(3, elements(1:4, tet)), [], 4)', tet);
axis equal
axis off
set(gcf, 'color', 'w');
% scatter3(x(s+1:maxNodes)', y(maxVertex+1:maxNodes)', z(maxVertex+1:maxNodes)');
% text(x+displacement, y+displacement, z+ displacement, c);
% patch('faces', edges, 'vertices', nodes','FaceVertexCData',edge_colors, 'FaceColor', 'flat');
hold off

end