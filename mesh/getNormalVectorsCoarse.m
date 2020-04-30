function [withNormals, bdryIdx, intIdx, bdry]= getNormalVectorsCoarse(mesh)

elements = mesh.Elements;
nodes = mesh.Nodes;
nodes = nodes - mean(nodes, 2);

temp = size(nodes);
maxNodes = temp(1, 2);
temp = size(elements);
numTets = temp(1,2);
maxVertex = max(elements(1:4, :), [], 'all');

% List of faces in counterclockwise order.
faceIndex = [
    3, 2, 1;
    2, 3, 4;
    1, 2, 4;
    4, 3, 1;
    ]';
faceIndex = reshape(faceIndex, [], 1);
faces = elements(faceIndex, 1:numTets);
flattened = reshape(faces, [], 1);
faces = reshape(flattened, 3, []);
sortedFaces = sort(faces)';
[C, ia, ic] = unique(sortedFaces, 'stable', 'rows');
tally = accumarray(ic, 1);
uniqueIndexes = find(tally - 2);
borderFaces = faces(:, ia(uniqueIndexes, 1))';

normals = cross(nodes(:, borderFaces(:, 1))-nodes(:, borderFaces(:, 2)), nodes(:, borderFaces(:, 1)) - nodes(:, borderFaces(:, 3)))';
normals = normalize(normals);
withNormals = [borderFaces(:, 1), normals, borderFaces(:, 2), normals, borderFaces(:, 3), normals]';
withNormals = reshape(withNormals, 4, [])';
withNormals = sortrows(withNormals);
indexes = withNormals(:, 1);


toKroneckerDelta1 = repmat(indexes, 1, maxVertex);
temp = size(indexes);
toKroneckerDelta2 = repmat(1:maxVertex, temp(1, 1), 1);
deltas = (toKroneckerDelta1 == toKroneckerDelta2);
deltas2 = deltas.*withNormals(:, 2);
deltas3 = deltas.*withNormals(:, 3);
deltas4 = deltas.*withNormals(:, 4);

[~, firstFaceIdx] = max(abs(deltas2) + abs(deltas2) + abs(deltas2));
firstFaceIdxes = sub2ind(size(deltas2), firstFaceIdx, 1:maxVertex);

withNormals = [deltas2(firstFaceIdxes); deltas3(firstFaceIdxes); deltas4(firstFaceIdxes)];
norms = vecnorm(withNormals);
bdryIdx = find(norms)';
intIdx = find(~norms)';
withNormals = (withNormals./norms)';

idxToBdryIdx = cumsum(norms > 0);

bdry = triangulation(reshape(idxToBdryIdx(borderFaces), [], 3), nodes(:, bdryIdx)');
end