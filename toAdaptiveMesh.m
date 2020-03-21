elements = model.Mesh.Elements;
nodes = model.Mesh.Nodes;
nodes = nodes - mean(nodes, 2);

[lambda, norms] = getLambda(model, mesh, qOcta);
tetCosts = mean(norms(elements))';

[~, idx] = sort(tetCosts(:, :), 'descend');
numTets = size(idx);
numTets = numTets(1, 1);
topTenPercent = idx(1:floor(numTets/8));
vertexIndex = [1, 2, 3, 4];
vertices = elements(vertexIndex, topTenPercent);
[vertices, ~, ~] = unique(vertices);

edgeIndex = [
    1, 2;
    1, 3;
    1, 4;
    2, 3;
    2, 4;
    3, 4;
    ]';
edgeIndex = reshape(edgeIndex, [], 1);
edges = elements(edgeIndex, topTenPercent);
flattened = reshape(edges, [], 1);
edges = reshape(flattened, 2, []);
sortedEdges = sort(edges)';
edges = unique(sortedEdges, 'stable', 'rows');

faceIndex = [
    3, 2, 1;
    2, 3, 4;
    1, 2, 4;
    4, 3, 1;
    ]';
faceIndex = reshape(faceIndex, [], 1);
faces = elements(faceIndex, topTenPercent);
flattened = reshape(faces, [], 1);
faces = reshape(flattened, 3, []);
sortedFaces = sort(faces)';
[faces, ia, ic] = unique(sortedFaces, 'stable', 'rows');


tetVertex = [elements; 1:numTets];
vertexIndex = [
    11, 1;
    11, 2;
    11, 3;
    11, 4;
]';

% tetVertex = tetVertex(:, tetVertex);
tetVertex = reshape(tetVertex(vertexIndex, :), 2, []);

vertexMat = zeros(numTets, maxNodes);
vertexMat(sub2ind(size(vertexMat), tetVertex(1, :), tetVertex(2, :))) = 1;
vertexMat = vertexMat(:, vertices);

numFaces = size(faces);
numFaces = numFaces(1, 1);
startFace = maxVertex;

numEdges = size(edges);
numEdges = numEdges(1, 1);
startEdge = startFace + numFaces;

faceToPoint = containers.Map(strtrim(cellstr(num2str(faces, '%d %d %d'))),  (startFace+1):(startFace+numFaces));
edgeToPoint = containers.Map(strtrim(cellstr(num2str(edges, '%d %d'))),  (startEdge + 1):(startEdge+numEdges));

tetsToAddMiddlePoint = zeros(numTets, 1);
numTetMiddlePoints = 0;
originalTets = zeros(numTets);
for tet=1:numTets
    b = 0;
    tetVertices = sort(elements(1:4, tet));
    tetFaceIndex = [
        1, 2, 3;
        1, 3, 4;
        1, 2, 4;
        2, 3, 4;
        ];
    tetFaces = tetVertices(tetFaceIndex);
    for faceInd=1:4
        face = tetFaces(faceInd, :);
        if isKey(faceToPoint, num2str(face, '%d %d %d'))
            b = 1;
            break;
        end
    end
    if b        
        numTetMiddlePoints = numTetMiddlePoints + 1;
        tetsToAddMiddlePoint(numTetMiddlePoints)= tet;
    else
        originalTets(tet - numTetMiddlePoints) = tet;
    end
end
startTet = numEdges + startEdge;
tetsToAddMiddlePoint = tetsToAddMiddlePoint(1:numTetMiddlePoints);
originalTets = originalTets(1:numTets - numTetMiddlePoints);
tetToPoint = containers.Map(tetsToAddMiddlePoint, (startTet+1):(startTet + numTetMiddlePoints));


newElements = zeros(12*numTets, 4);
counter = 0;
for tet=1:numTets
    if isKey(tetToPoint, tet)
        tetVertices = sort(elements(1:4, tet));
        tetFaceIndex = [
            1, 2, 3;
            1, 3, 4;
            1, 2, 4;
            2, 3, 4;
            ];
        tetFaces = tetVertices(tetFaceIndex);
        tetVert = tetToPoint(tet);
        for faceInd=1:4
            faceEdgeIndex = [
                1, 2;
                1, 3;
                2, 3;
                ];
            face = tetFaces(faceInd, :);
            faceEdges = face(faceEdgeIndex);
            if isKey(faceToPoint, num2str(face, '%d %d %d'))
                faceVert = faceToPoint(num2str(face, '%d %d %d'));
                for edgeInd=1:3
                    edge = faceEdges(edgeInd, :);
                    edgeVert = edgeToPoint(num2str(edge, '%d %d'));
                    counter = counter + 1;
                    newElements(counter, :) = [tetVert, edgeVert, faceVert, edge(1, 1)]; 
                    counter = counter + 1;
                    newElements(counter, :) = [tetVert, edgeVert, faceVert, edge(1, 2)];
                end
            else
                counter = counter + 1;
                newElements(counter, :) = [tetVert, face(1), face(2), face(3)];
            end
        end
    else
        vertexEdgeIndex = [
            1, 2;
            1, 3;
            2, 3;
            1, 4;
            2, 4;
            3, 4;
            ];
        for edgeInd=1:6
            key = num2str(elements(vertexEdgeIndex(edgeInd, :), tet), '%d %d');
            if isKey(edgeToPoint, key)
                counter = counter + 1
                break
            end 
        end
    end
end

newElements = newElements(1:counter, :);
newElements = sort(newElements, 2);

edgeIndex = [
    1, 2;
    2, 3;
    1, 3;
    1, 4;
    2, 4;
    3, 4;
    ]';


newEdges = unique(reshape(newElements(:, edgeIndex)', 2, [])', 'stable', 'rows');
numNewEdges = size(newEdges);
numNewEdges = numNewEdges(1, 1);
startQuadratic = startTet + numTetMiddlePoints;
edgeToMidpoint = containers.Map(strtrim(cellstr(num2str(newEdges, '%d %d'))), (startQuadratic+1):(startQuadratic+numNewEdges));

newElements = [newElements'; zeros(6, size(newElements, 1))];
for tet=1:size(newElements, 2)
    edges = reshape(newElements(edgeIndex, tet), 2, [])';
    for i=1:6
        newElements(4+i, tet) = edgeToMidpoint(num2str(edges(i, :), '%d %d'));
    end
end

numNewNodes = startQuadratic+numNewEdges;

newNodes = zeros(numNewNodes, 3);
newNodes(1:maxVertex, :) = nodes(:, 1:maxVertex)';
for key=faceToPoint.keys
    verts = str2num(cell2mat(key));
    newNodes(faceToPoint(key{1,1}), :) = mean(newNodes(verts, :), 1);
end

for key=edgeToPoint.keys
    verts = str2num(cell2mat(key));
    newNodes(edgeToPoint(key{1,1}), :) = mean(newNodes(verts, :), 1);
end

for key=tetToPoint.keys
    verts = elements(1:4, key{1,1});
    newNodes(tetToPoint(key{1,1}), :) = mean(newNodes(verts, :), 1);
end

for key=edgeToMidpoint.keys
    verts = str2num(cell2mat(key));
    newNodes(edgeToMidpoint(key{1,1}), :) = mean(newNodes(verts, :), 1);
end

newNodes = newNodes';
newElements = reorientTets(newElements, newNodes);

visualizeMesh(newNodes, newElements);

newModel = createpde('thermal', 'transient');
geometryFromMesh(newModel, newNodes, newElements(1:4, :));

visualizeMesh(newNodes, newElements);

function [elements]=reorientTets(elements, nodes)
    numTets = size(elements);
    numTets = numTets(1, 2);
    
    switchOrientation = [
        2, 1, 3, 4, 5, 7, 6, 9, 8, 10
        ];
    
    for tet = 1:numTets
        vec1 = nodes(:, elements(2, tet)) - nodes(:, elements(1, tet));
        vec2 = nodes(:, elements(3, tet)) - nodes(:, elements(1, tet));
        normal = cross(vec1, vec2);
        signOfDot = dot(normal, nodes(:, elements(4, tet)) - nodes(:, elements(1, tet)));
%         visualizeMesh(nodes, elements(:, tet));
        if signOfDot > 0
            continue
        else
            elements(:, tet) = elements(switchOrientation, tet);
        end
    end
end

