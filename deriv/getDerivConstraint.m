mesh.nv = size(model.Mesh.Nodes, 2);

dg = getDerivOfG(mesh, qOcta);
d2g = getHessianOfG(mesh);
IHh = buildIHh(model.Mesh);
qFine = (IHh*qOcta')';
lambdaH = sparse(getLambdaH(mesh, qFine, dg));

numBdry = size(mesh.bdryIdx, 1);
numInt = size(mesh.intIdx, 1);

meshEdge = struct;
meshEdge.nv = mesh.nv;
meshEdge.bdryIdx = mesh.bdryIdxFine(numBdry+1:end, :);
meshEdge.intIdx = mesh.intIdxFine(numInt+1:end, :);
meshEdge.bdryNormals = mesh.bdryNormalsFine(numBdry+1:end, :);
dgEdge = getDerivOfG(meshEdge, qFine);

dotProd = lambdaH*d2g';
dotProd = reshape(dotProd, mesh.nv*9, []);

d2f = getHessianOfF(mesh, qOcta);

numConstraintsCoarse = size(dg, 2);
numConstraintsFine = size(dgEdge, 2);
numNodes = mesh.nv * 9;
numTotal = numNodes + numConstraintsCoarse + numConstraintsFine;
constMat = [d2f+dotProd, dg, dgEdge; dg', sparse(numConstraintsCoarse, numConstraintsCoarse), sparse(numConstraintsCoarse,numConstraintsFine)];
numNewConsts = size(constMat, 1);

B = [speye(numNodes), sparse(numNodes, numConstraintsFine + numConstraintsCoarse)];
x0 = getDerivOfF(mesh, qFine);
C = [B'*sparse(x0); sparse(numNewConsts, 1)];
M = [2 * (B'*B); constMat];

xComb = C \ M;
numIntFine = size(meshEdge.intIdx, 1);
numBdryFine = size(meshEdge.bdryIdx, 1);

dLambda = xComb((numTotal - numConstraintsFine + 1):numTotal);
dLambdaBdry = reshape(dLambda(1:(numBdryFine*8))', numBdryFine, []);
dLambdaInt = reshape(dLambda((numBdryFine*8 + 1):(numBdryFine*8+numIntFine*15))', numIntFine, []);

edgeCosts = zeros(mesh.nv, 1);
edgeCosts(meshEdge.bdryIdx) = vecnorm(dLambdaBdry, 2, 2);
edgeCosts(meshEdge.intIdx) = vecnorm(dLambdaInt, 2, 2);