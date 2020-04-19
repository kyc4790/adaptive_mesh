function IHh=buildIHh(mesh)
elements = mesh.Elements;
nodes = mesh.Nodes;
nodes = nodes - mean(nodes, 2);

temp = size(nodes);
maxNodes = temp(1, 2);
maxVertex = max(elements(1:4, :), [], 'all');

edgeEndpointIndex1 = [1; 2; 1; 1; 2; 3]';
edgeEndpointIndex2 = [2; 3; 3; 4; 4; 4]';
edgeEndpoint1 = elements(edgeEndpointIndex1, :);
edgeEndpoint2 = elements(edgeEndpointIndex2, :);
IHh = zeros(maxNodes, maxVertex);

sz = [maxNodes, maxVertex];
oneInd = sub2ind(sz, 1:maxVertex, 1:maxVertex);
IHh(oneInd) = 1;

halfInd = sub2ind(sz, elements(5:10, :), edgeEndpoint1);
IHh(halfInd) = 0.5;

halfInd = sub2ind(sz, elements(5:10, :), edgeEndpoint2);
IHh(halfInd) = 0.5;

end