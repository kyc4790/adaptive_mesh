function MBOmesh = buildMesh(mesh)
    MBOmesh = struct;
    nodes = mesh.Nodes;
    elements = mesh.Elements;
    nodes = nodes - mean(nodes, 2);
    maxVertex = max(elements(1:4, :), [], 'all');
    temp = size(nodes);
    maxNodes = temp(1, 2);
    MBOmesh.nv = maxNodes;
    [withNormals, bdryIdx, intIdx, bdry] = getNormalVectorsFine(mesh);
    MBOmesh.bdryIdx = bdryIdx;
    MBOmesh.intIdx = intIdx;
    MBOmesh.bdryNormals = withNormals(bdryIdx,:);
    MBOmesh.bdry = bdry;
    [withNormals, bdryIdx, intIdx, bdry] = getNormalVectorsFine(mesh);
    MBOmesh.bdryIdxFine = bdryIdx;
    MBOmesh.intIdxFine = intIdx;
    MBOmesh.bdryNormalsFine = withNormals(bdryIdx,:);
    MBOmesh.bdryFine = bdry;
    
    model = createpde('thermal', 'transient');
    geometryFromMesh(model, mesh.Nodes, mesh.Elements);
    thermalProperties(model,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);
    
    FEM = assembleFEMatrices(model);
    MBOmesh.L = FEM.K; %(1:maxVertex, 1:maxVertex);
    MBOmesh.M = FEM.M; %(1:maxVertex, 1:maxVertex);
    MBOmesh.fineL = FEM.K;
    
    [eigvectors, eigvalues] = eigs(MBOmesh.L, 2, 'sm');
    MBOmesh.lambda1L = eigvalues(2, 2);
    MBOmesh.eigvalues = eigvalues;
    MBOmesh.eigvectors = eigvectors;
    
    MBOmesh.tetra = triangulation(mesh.Elements(1:4, :)', nodes');
end