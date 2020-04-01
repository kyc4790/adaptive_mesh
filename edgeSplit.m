function [newModel, newQ]=edgeSplit(model, costs, qOcta)
%% Setup
    elements = model.Mesh.Elements;
    nodes = model.Mesh.Nodes;
    nodes = nodes - mean(nodes, 2);

    maxVertex = max(elements(1:4, :),[], 'all');
    maxNodes = max(elements,[], 'all');
    [~, idx] = sort(costs(maxVertex+1:maxNodes));

    numTets = size(elements);
    numTets = numTets(1, 2);

    %%  Finds set of edges to be split
    tetSet = zeros(1, numTets);
    index = 1;
    numEdgesPicked = 0;
    RATIO = 0.05;
    edgesToPick = ceil(RATIO * (maxNodes - maxVertex));
    pickedEdges = zeros(edgesToPick, 1);
    oldToNew = 1:maxNodes;
    while numEdgesPicked < edgesToPick && index <= maxNodes - maxVertex
        edge = idx(index) + maxVertex;
        [overlap, tetSet] = hasOverlap(tetSet, edge, elements);
        if overlap==1
        else
            numEdgesPicked = numEdgesPicked + 1;
            pickedEdges(numEdgesPicked) = edge;
            temp = oldToNew(edge);
            oldToNew(edge) = oldToNew(maxVertex + numEdgesPicked);
            oldToNew(maxVertex + numEdgesPicked) = temp;
        end
        index = index + 1;
    end

    %% Constructs new mesh + new q
    newTets = sum(tetSet, 'all');

    newQ = qOcta(:, oldToNew(1:(maxVertex + numEdgesPicked)));
    newNodes = [nodes(:, oldToNew), zeros(3, numEdgesPicked * 2)];
    newElements = [oldToNew(elements(:, find(1-tetSet)))'; zeros(10, 2 * newTets)']';

    index = numTets - newTets + 1;
    for i=1:numEdgesPicked
        % not quite enough edges
        if pickedEdges(i) == 0
            break
        end

        edge = pickedEdges(i);
        newEdge = oldToNew(edge);
        [idx, tets] = find(elements==edge);
        numNewTets = size(tets);
        numNewTets = numNewTets(1, 1);

        for j=1:numNewTets
            tet = tets(j);
            endpoints = [0, 0];
            switch idx(j)
                case 5 
                    endpoints = [1, 2];
                case 6 
                    endpoints = [2, 3];
                case 7 
                    endpoints = [1, 3];
                case 8 
                    endpoints = [1, 4];
                case 9 
                    endpoints = [2, 4];
                case 10 
                    endpoints = [3, 4];
            end
            if elements(endpoints(1), tet) > elements(endpoints(2), tet)
                endpoints = flip(endpoints);
            end
            

            newElements(:, index) = oldToNew(elements(:, tet));
            newElements(:, index+1) = oldToNew(elements(:, tet));
            newElements(endpoints(1), index) = newEdge;
            newElements(endpoints(2), index+1) = newEdge;
            newElements(idx(j), index) = 2*i + maxNodes - 1;
            newElements(idx(j), index + 1) = 2*i + maxNodes;
            
            endpoints = elements(endpoints, tet);
            newNodes(:, 2*i + maxNodes - 1) = mean([nodes(:, edge), nodes(:, endpoints(1))], 2);
            newNodes(:, 2*i + maxNodes) = mean([nodes(:, edge), nodes(:, endpoints(2))], 2);
            index = index + 2;
        end

    end
    newModel = createpde('thermal', 'transient');
    visualizeMesh(newNodes, newElements);
    geometryFromMesh(newModel, newNodes, newElements);
end

function [overlap, tetSet]=hasOverlap(tetSet, edge, elements)
    %% Checks for overlaps
    % checks that no tets in tetSet have an edge
    % if none, adds all tets bordering the edge to tetSet
    numTets = sum(tetSet, 'all');
    tets = find(tetSet);
    
    subMat = elements(:, tets);
    
    numNotEdge = size(find(subMat==edge));
    numNotEdge = numNotEdge(1, 1);
    if numNotEdge == 0
        overlap = 0;
        [~, idx] = find(elements==edge);
        tetSet(idx) = 1;            
    else
        overlap = 1;
    end   
end