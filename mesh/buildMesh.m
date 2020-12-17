function MBOmesh = buildMesh(mesh, isFine)
    if nargin < 2
        isFine = true;
    end
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
    if(isFine)
        [withNormals, bdryIdx, intIdx, bdry] = getNormalVectorsFine(mesh);
        MBOmesh.bdryIdxFine = bdryIdx;
        MBOmesh.intIdxFine = intIdx;
        
        model = createpde('thermal', 'transient');
        geometryFromMesh(model, mesh.Nodes, mesh.Elements);
        thermalProperties(model,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);
        
        FEM = assembleFEMatrices(model);
        MBOmesh.fineL = FEM.K;
        MBOmesh.fineM = FEM.M;
    end
    MBOmesh.bdryNormalsFine = withNormals(bdryIdx,:);
    MBOmesh.bdryFine = bdry;
    
    modelCoarse = createpde('thermal', 'transient');
    geometryFromMesh(modelCoarse, mesh.Nodes(:, 1:maxVertex), mesh.Elements(1:4, :));
    thermalProperties(modelCoarse,'ThermalConductivity',0.08, 'MassDensity', 1, 'SpecificHeat', 1);
    
    FEMCoarse = assembleFEMatrices(modelCoarse);
    MBOmesh.L = FEMCoarse.K; %(1:maxVertex, 1:maxVertex);
    MBOmesh.M = FEMCoarse.M; %(1:maxVertex, 1:maxVertex);

    
    [eigvectors, eigvalues] = eigs(MBOmesh.L, 2, 'sm');
    MBOmesh.lambda1L = eigvalues(2, 2);
    MBOmesh.eigvalues = eigvalues;
    MBOmesh.eigvectors = eigvectors;
    
    warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    MBOmesh.tetra = triangulation(mesh.Elements(1:4, :)', nodes');
    warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');

end