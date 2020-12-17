% Checks that projBoundary functions correctly by attempting to project a
% vector onto the boundary.
%

alpha = 0.01;

b = [0 0 0 sqrt(7/12) 0 0 0]';
numBdry = size(weightedMesh.bdryNormals, 1);
B = repmat(b, 1, numBdry);

qProj = projBoundary(alpha, weightedMesh.bdryNormals, ones(9, numBdry));

mat = OctaAlignMat(weightedMesh.bdryNormals);
W = multiprod(multitransp(mat), sparse(2:8, 1:7, ones(1, 7), 9, 7));
res = multiprod(multitransp(W), qProj);
q_inds = repelem(1:numBdry, 7, 1);
const_inds = repelem(1:7, numBdry, 1)';
inds = sub2ind(size(res), const_inds, q_inds, q_inds);
vals = res(inds) - B;
assert(max(max(abs(vals))) <= alpha + 1e6, 'not within threshold')