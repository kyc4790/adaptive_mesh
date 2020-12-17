function MBOmesh = buildWeightedMesh(mesh)
    MBOmesh = struct;
    nodes = mesh.Nodes;
    elements = mesh.Elements;
    nodes = nodes - mean(nodes, 2);
    maxVertex = max(elements(1:4, :), [], 'all');
    temp = size(nodes);
    maxNodes = temp(1, 2);
    MBOmesh.nv = maxVertex;
    [withNormals, bdryIdx, intIdx, bdry] = getNormalVectorsCoarse(mesh);
    MBOmesh.bdryIdx = bdryIdx;
    MBOmesh.intIdx = intIdx;
    MBOmesh.bdryNormals = withNormals(bdryIdx,:);
    MBOmesh.bdry = bdry;
    
    MBOmesh.weights = zeros(maxVertex, 1);
    MBOmesh.weights(bdryIdx, 1) = 1;
    
    edgeIndexes = [1, 2; 2, 3; 1, 3; 1, 4; 2, 4; 3, 4]';
    edges = reshape(elements(edgeIndexes, :), 2, []);
    edges = unique(sort(edges)', 'rows');
    adjMat = sparse(edges(:, 1), edges(:, 2), vecnorm(nodes(:, edges(:, 1)) - nodes(:, edges(:, 2))), maxVertex, maxVertex);
    adjMat = adjMat + adjMat';
    
    dists = graphallshortestpaths(adjMat);    
    
    bdryDists = dists(bdryIdx, :);
    MBOmesh.weights = min(bdryDists) + 0.1;
    
    [Lw, L, M, ~] = ComputeWeightedLaplacian(nodes(:, 1:maxVertex)', elements(1:4, :)', MBOmesh.weights);
%     model = createpde('thermal', 'transient');
%     modelCoarse = createpde('thermal', 'transient');
%     geometryFromMesh(model, mesh.Nodes, mesh.Elements);
%     geometryFromMesh(modelCoarse, mesh.Nodes(:, 1:maxVertex), mesh.Elements(1:4, :));
%     thermalProperties(model,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);
%     thermalProperties(modelCoarse,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);
    
    MBOmesh.L = L;
    MBOmesh.M = M;
    
    [eigvectors, eigvalues] = eigs(MBOmesh.L, 2, 'sm');
    MBOmesh.lambda1L = eigvalues(2, 2);
    MBOmesh.eigvalues = eigvalues;
    MBOmesh.eigvectors = eigvectors;
    
    MBOmesh.L = Lw;
    
    warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    MBOmesh.tetra = triangulation(mesh.Elements(1:4, :)', nodes');
    warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');

end