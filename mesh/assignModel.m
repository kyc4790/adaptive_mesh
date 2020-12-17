function newModel=assignModel(mesh)
    % Builds a quadratic FEModel
    % @param mesh   an FEMesh
    %

    elements = mesh.Elements;
    maxVertex = max(max(elements(1:4, :)));
    nodes = mesh.Nodes(:, 1:maxVertex);
    numNodes = size(nodes, 2);

    % get all the edges
    edgeIndexes = [1, 2; 2, 3; 1, 3; 1, 4; 2, 4; 3, 4]';
    edges = reshape(elements(edgeIndexes, :), 2, []);
    edges = unique(sort(edges)', 'rows');
   
    % get locations of everything
    all_locs = reshape(nodes(:, edges'), 6, []);
    edgeNodes = (all_locs(1:3, :) + all_locs(4:6, :)) / 2;
    
    % assign each edge a value
    edgeLocs = sparse(numNodes, numNodes);
    edgeLocs(sub2ind(size(edgeLocs),  edges(:, 1), edges(:, 2))) = ((1:size(edges, 1)) + numNodes); 
    edgeLocs = max(edgeLocs, edgeLocs');
    
    % get tetrahedral values
    newElements = [elements; zeros(6, size(elements, 2))];
    newElements(5, :) = edgeLocs(sub2ind(size(edgeLocs), elements(1, :), elements(2,:)));
    newElements(6, :) = edgeLocs(sub2ind(size(edgeLocs), elements(2, :), elements(3,:)));
    newElements(7, :) = edgeLocs(sub2ind(size(edgeLocs), elements(1, :), elements(3,:)));
    newElements(8, :) = edgeLocs(sub2ind(size(edgeLocs), elements(1, :), elements(4,:)));
    newElements(9, :) = edgeLocs(sub2ind(size(edgeLocs), elements(2, :), elements(4,:)));
    newElements(10, :) = edgeLocs(sub2ind(size(edgeLocs), elements(3, :), elements(4,:)));
    
    newNodes = [nodes, edgeNodes];
    
    newModel = createpde('thermal', 'transient');
    geometryFromMesh(newModel, newNodes, newElements);

end