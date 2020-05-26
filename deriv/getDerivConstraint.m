function [edgeCosts, qFine, dLambda]=getDerivConstraint(model, mesh, qOcta)
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
% 
%     dotProd = lambdaH*d2g';
%     dotProd = reshape(dotProd, mesh.newnv*9, []);
% 
%     d2f = getHessianOfF(mesh, qOcta);
% 
%     numConstraintsCoarse = size(dg, 2);
%     numConstraintsFine = size(dgEdge, 2);
%     numNodes = mesh.newnv * 9;
%     numTotal = numNodes + numConstraintsCoarse + numConstraintsFine;
%     constMat = [d2f+dotProd, dg, dgEdge; 
%         dg', sparse(numConstraintsCoarse, numConstraintsCoarse), sparse(numConstraintsCoarse,numConstraintsFine)];
%     numNewConsts = size(constMat, 1);
% 
%     B = [speye(numNodes), sparse(numNodes, numConstraintsFine + numConstraintsCoarse)];
%     x0 = getDerivOfF(mesh, qFine);
%     C = [B'*sparse(x0); sparse(numNewConsts, 1)];
%     M = [(B'*B), constMat'; constMat, sparse(numNewConsts, numNewConsts)];
%     
%     xComb = M \ C;
%     
%     A = B; b = x0;
%     C = constMat;
%     cvx_begin quiet
%         cvx_solver mosek
%         cvx_precision high
%         variable x(numTotal)
%         minimize( norm( A * x - b, 2 ) )
%         subject to
%             C * x == 0
%     cvx_end
% 
% %     cvx_begin
% %         cvx_precision best
% %         cvx_solver mosek
% %         variable xComb(numTotal)
% %         minimize( norm( B * xComb - x0, 2 ) )
% %         subject to
% %             constMat * xComb == 0
% %     cvx_end
%     
%     xComb = x;
    
    numIntFine = size(meshEdge.intIdx, 1);
    numBdryFine = size(meshEdge.bdryIdx, 1);

    dLambda = getG(meshEdge, qFine);
%     dLambda = xComb((numTotal - numConstraintsFine + 1):numTotal);
%     
%     checkEnvelopeTheorem(meshEdge, qFine, dLambda);
    
    dLambdaBdry = reshape(dLambda(1:(numBdryFine*8))', numBdryFine, []);
    dLambdaInt = reshape(dLambda((numBdryFine*8 + 1):(numBdryFine*8+numIntFine*15))', numIntFine, []);

    edgeCosts = zeros(mesh.newnv, 1);
    edgeCosts(meshEdge.bdryIdx) = vecnorm(dLambdaBdry, 2, 2);
    edgeCosts(meshEdge.intIdx) = vecnorm(dLambdaInt, 2, 2);
end
