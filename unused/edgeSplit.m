function [newModel, newQ, newNodes, newElements]=edgeSplit(model, costs, qOcta)
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
    newNodes = [nodes(:, oldToNew), zeros(3, numEdgesPicked * 3 + newTets)];
    newElements = [oldToNew(elements(:, find(1-tetSet)))'; zeros(10, 2 * newTets)']';

    index = numTets - newTets + 1;
    edgeIndex = maxNodes;
    for i=1:numEdgesPicked
        if pickedEdges(i) == 0
            break
        end

        edge = pickedEdges(i);
        newEdge = oldToNew(edge);
        [idx, tets] = find(elements==edge);
        numNewTets = size(tets);
        numNewTets = numNewTets(1, 1);
        vertices = unique(elements(1:4, tets)); % unique mapping!
        numNewEdges = size(vertices, 1);
        
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

            for endpoint=1:2
                newElements(:, index) = oldToNew(elements(:, tet));
                newElements(endpoints(endpoint), index) = newEdge;
                correspondingVertex = 1:4;
                correspondingVertex = find(correspondingVertex - endpoints(endpoint));                
                correspondingEdge = zeros(1, 3);
                switch(endpoints(endpoint))
                    case 1
                        correspondingEdge=[5, 7, 8];
                    case 2
                        correspondingEdge=[5, 6, 9];
                    case 3
                        correspondingEdge=[6, 7, 10];
                        correspondingVertex=[2, 1, 4];
                    case 4
                        correspondingEdge=[8, 9, 10];
                end
                for elem=1:3
                    edge_elem = find(vertices == elements(correspondingVertex(elem), tet));
                    newElements(correspondingEdge(elem), index) = edgeIndex + edge_elem;
                end
                index = index + 1;
            end            
        end
        
        for j=1:numNewEdges
            edgeIndex = edgeIndex + 1;
            newNodes(:, edgeIndex) = mean([nodes(:, vertices(j)), nodes(:, edge)], 2);
        end
    end
    newNodes=newNodes(:, 1:edgeIndex);
    newModel = createpde('thermal', 'transient');
    visualizeMesh(newNodes, newElements);
    try
        geometryFromMesh(newModel, newNodes, newElements);
    catch ME1
        for iter=1:(numTets + newTets)
            newModel = createpde('thermal', 'transient');
            try
                geometryFromMesh(newModel, newNodes, newElements(:, iter));
            catch ME2
                disp(iter);
                disp(newElements(:, iter));
            end
        end
        disp(ME1);
    end
end

function [overlap, tetSet]=hasOverlap(tetSet, edge, elements)
    %% Checks for overlaps
    % checks that no tets in tetSet have an edge
    % if none, adds all tets bordering the edge to tetSet
    tets = find(tetSet);
    
    subMat = elements(:, tets);
    
    numEdge = size(find(subMat==edge));
    numEdge = numEdge(1, 1);
    if numEdge == 0
        overlap = 0;
        [~, idx] = find(elements==edge);
        tetSet(idx) = 1;            
    else
        overlap = 1;
    end   
end