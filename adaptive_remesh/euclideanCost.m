function [costs, qFine] = euclideanCost(q, model)
% Computes the cost of edges based on Euclidean distance in R^9
%
% @param q          the original q to be processed
% @param model      the model associated with the tet mesh.
%
% @returns costs    an (num_vertices + num_edges) by 1 vector of costs
% @returns qFine    the interpolated value of q based on the quadratic
% model
    % Get unique edges
   edgeIndexes = [5, 1, 2; 6, 2, 3;7, 1, 3;8, 1, 4;9, 2, 4;10, 3, 4]';
   edges = reshape(model.Mesh.Elements(edgeIndexes, :), 3, []);
   edges = unique(sort(edges)', 'rows');
   edges = [edges(:, 3), edges(:, 1:2)];
   edges = sortrows(edges);
   % Get euclidean distance
   val1 = q(:, edges(:, 2));
   val2 = q(:, edges(:, 3));
   
   loc1 = model.Mesh.Nodes(:, edges(:, 2));
   loc2 = model.Mesh.Nodes(:, edges(:, 3));
   
   norms = vecnorm(val1 - val2);
   dists = vecnorm(loc1 - loc2);
   
   costs = [zeros(max(max(model.Mesh.Elements(1:4, :))), 1); (norms)'];
   
   IHh = buildIHh(model.Mesh);
   qFine = (IHh*q')';
   
end