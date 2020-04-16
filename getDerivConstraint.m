dg = getDerivOfG(mesh, qOcta);
d2g = getHessianOfG(mesh);
lambdaH = sparse(getLambdaH(mesh, qOcta, dg));

dotProd = lambdaH*d2g';
dotProd = reshape(dotProd, mesh.nv*9, []);

d2f = sparse(9*mesh.nv, 9*mesh.nv);

inds = ((1:mesh.nv)-1)*9;
for i=1:9
     d2f(inds+i, inds+i) = mesh.L;
end

numConstraintsCoarse = size(dg, 2);
constMat = [d2f+dotProd, dg; dg', sparse(numConstraintsCoarse, numConstraintsCoarse)];