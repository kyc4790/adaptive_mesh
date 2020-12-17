function [mesh]=getModel(gmsh)
    % Builds a FEMesh from a gmsh generated MATLAB file
    % gmsh the object loaded from the gmsh generated MATLAB file
    
    tets = gmsh.TETS(:, 1:4)';
    nodes = gmsh.POS';
    mesh.Nodes = nodes;
    mesh.Elements = tets;
    
%     % get the edges
%     edgeIndexes = [1, 2; 2, 3; 1, 3; 1, 4; 2, 4; 3, 4]';
%     edges = reshape(tets(edgeIndexes, :), 2, []);
%     edges = unique(sort(edges)', 'rows');
    
%     % Assign the edges a number
%     edge_nums = [(gmsh.nbNod+(1:size(edges, 1)))', edges];
%     nodes = [nodes, zeros(3, size(edges, 2))];
%     
%     nodes(1, edge_nums(:, 1)) = mean(reshape(nodes(1, edge_nums(:, 2:3)), [], 2), 2);
%     nodes(2, edge_nums(:, 1)) = mean(reshape(nodes(1, edge_nums(:, 2:3)), [], 2), 2);
%     nodes(2, edge_nums(:, 1)) = mean(reshape(nodes(1, edge_nums(:, 2:3)), [], 2), 2);

    
end