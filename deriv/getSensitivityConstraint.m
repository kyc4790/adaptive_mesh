function [edgeCosts, qFine, dLambda]=getSensitivityConstraint(model, mesh, qOcta)
    mesh.newnv = size(model.Mesh.Nodes, 2);
    dg = getDerivOfG(mesh, qOcta);
    d2g = getHessianOfG(mesh);

    IHh = buildIHh(model.Mesh);
    qFine = (IHh*qOcta')';
    lambdaH = sparse(getLambdaH(mesh, qFine, dg));
    
    numBdry = size(mesh.bdryIdx, 1);
    numInt = size(mesh.intIdx, 1);
    
    meshEdge = struct;
    meshEdge.newnv = mesh.newnv;
    meshEdge.bdryIdx = mesh.bdryIdxFine(numBdry+1:end, :);
    meshEdge.intIdx = mesh.intIdxFine(numInt+1:end, :);
    meshEdge.bdryNormals = mesh.bdryNormalsFine(numBdry+1:end, :);
    dgEdge = getDerivOfG(meshEdge, qFine);
    
    dotProd = lambdaH*d2g';
    dotProd = reshape(dotProd, mesh.newnv*9, []);
    
    df = getDerivOfF(mesh, qFine);
    d2f = getHessianOfF(mesh, qFine);
    
    xLength = 9 * mesh.newnv;
    numConstraintsCoarse = size(dg, 2);
    numConstraintsFine = size(dgEdge, 2);
    
    Cy = [d2f + dotProd, dg; dg', sparse(numConstraintsCoarse, numConstraintsCoarse)];
    Cz = [dgEdge; sparse(numConstraintsCoarse, numConstraintsFine)];
        
    A = Cy \ Cz;
    Jy = [df; sparse(numConstraintsCoarse, 1)];
    
    dLambda = -Jy'*A;
    checkEnvelopeTheorem(meshEdge, qFine, dLambda);
    
    dLambdaBdry = reshape(dLambda(1:(numBdryFine*8))', numBdryFine, []);
    dLambdaInt = reshape(dLambda((numBdryFine*8 + 1):(numBdryFine*8+numIntFine*15))', numIntFine, []);

    edgeCosts = zeros(mesh.newnv, 1);
    edgeCosts(meshEdge.bdryIdx) = vecnorm(dLambdaBdry, 2, 2);
    edgeCosts(meshEdge.intIdx) = vecnorm(dLambdaInt, 2, 2);
end