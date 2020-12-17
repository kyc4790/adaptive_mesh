function MBOmesh = adjustWeights(model, mesh, fun)
% Returns MBOmesh, a mesh that can be passed into MBO with weights in
% fun(mesh.weights)
% 
% @param model the triangulation
% @param mesh  the original weighted MBO mesh that contains originally
%              calculated weights
% @param fun   the function to be applied to the weights in mesh

    elements = model.Elements;
    nodes = model.Nodes;
    maxVertex = mesh.nv;
    MBOmesh.nv = maxVertex;
    MBOmesh.bdryIdx = mesh.bdryIdx;
    MBOmesh.intIdx = mesh.intIdx;
    MBOmesh.bdryNormals = mesh.bdryNormals;
    MBOmesh.bdry = mesh.bdry;
    
    MBOmesh.weights = fun(mesh.weights);
    
    MBOmesh.lambda1L = mesh.lambda1L;
    MBOmesh.eigvalues = mesh.eigvalues;
    MBOmesh.eigvectors = mesh.eigvectors;
    
    [Lw, ~, M, ~] = ComputeWeightedLaplacian(nodes(:, 1:maxVertex)', elements(1:4, :)', MBOmesh.weights);
    MBOmesh.L = Lw;
    MBOmesh.M = M;
    
    MBOmesh.tetra = mesh.tetra;
end