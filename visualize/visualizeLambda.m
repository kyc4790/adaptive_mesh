% IHh = buildIHh(model.Mesh);
% 
% A = [eye(maxNodes); IHh'];
% B = [2*mesh.fineL*IHh*qOcta'; (mesh.L * qOcta.')];
% 
% lambda = B\A;
% norms = vecnorm(lambda);

x = nodes(1, :)';
y = nodes(2, :)';
z = nodes(3, :)';
plotData = struct();
plotData.x = x;
plotData.y = y;
plotData.z = z;
plotData.colors = costs;
scatter3(x, y, z, 100, plotData.colors, 'filled');
hold on
axis equal
axis off
set(gcf, 'color', 'w');
hold off