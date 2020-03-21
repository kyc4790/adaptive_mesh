function [lambda, norms]=getLambda(model, mesh, qOcta)
IHh = buildIHh(model.Mesh);

maxNodes = size(mesh.fineL);
maxNodes = maxNodes(1);
A = [eye(maxNodes); IHh'];
B = [2*mesh.fineL*IHh*qOcta'; (mesh.L * qOcta.')];

lambda = B\A;
norms = vecnorm(lambda);
end