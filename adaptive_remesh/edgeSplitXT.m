function [newModel, newQ]=edgeSplitXT(model, costs, qOcta, mesh, RATIO)%
% Detects the edges to split and splits them
% Returns an interpolated value of q corresponding to the subdivided mesh
%
% @param model      the model associated with the mesh
% @param costs      an (num_vertices + num_edges) by 1 vector of the costs associated with each control point 
% @param qOcta      the original qOcta
% @param mesh       the MBOmesh associated with the mesh
% @param RATIO      the ratio of edges to divide (default 0.05)
%
% @returns newModel the model associated with subdivided mesh
% @returns newQ     the q associated with the subdivided mesh

    %% Setup
    elements = model.Mesh.Elements;
    nodes = model.Mesh.Nodes;
    nodes = nodes - mean(nodes, 2);

    maxVertex = max(elements(1:4, :),[], 'all');
    maxNodes = max(elements,[], 'all');
    
    edgeIndexes = [1, 2, 5; 2, 3, 6; 1, 3, 7; 1, 4, 8; 2, 4, 9; 3, 4, 10]';
    edges = reshape(elements(edgeIndexes, :), 3, []);
    edges = unique(sort(edges)', 'rows');
    edges = sortrows([edges(:, 3:3), edges(:, 1:2)]);
    % (first vertex, second vertex, edge index)
    bdryIdx = mesh.bdryIdx;
    isBdry = all(ismember(edges(:,2:3), bdryIdx), 2);
    
    [~, idx] = sort(costs(maxVertex+1:maxNodes), 'descend');

    numTets = size(elements);
    numTets = numTets(1, 2);

    %%  Finds set of edges to be split
    tetSet = zeros(1, numTets);
    index = 1;
    numEdgesPicked = 0;
    
    if(nargin < 5)
        RATIO = 0.05;
    end
    
    
    edgesToPick = ceil(RATIO * (maxNodes - maxVertex));
    pickedEdges = zeros(edgesToPick, 1);
    oldToNew = 1:maxNodes;
    newToOld = 1:maxNodes;
    while numEdgesPicked < edgesToPick && index <= maxNodes - maxVertex
        edge = idx(index) + maxVertex;
%         if(isBdry(idx(index), 1))
%             index = index + 1;
%             continue;
%         end
        [overlap, tetSet] = hasOverlap(tetSet, edge, elements);
        if overlap==1
        else
            numEdgesPicked = numEdgesPicked + 1;
            pickedEdges(numEdgesPicked) = edge;
            
            oldToNew(oldToNew==edge) = oldToNew(maxVertex + numEdgesPicked);
            oldToNew(maxVertex + numEdgesPicked) = edge;
            
            newToOld(newToOld==(maxVertex + numEdgesPicked)) = newToOld(edge);
            newToOld(edge) = maxVertex + numEdgesPicked;
        end
        index = index + 1;
    end

    %% Constructs new mesh + new q
    newTets = sum(tetSet, 'all');
    
    newQ = qOcta(:, oldToNew(1:(maxVertex + numEdgesPicked)));
    newNodes = nodes(:, oldToNew(1:(maxVertex + numEdgesPicked)));
    newElements = [newToOld(elements(1:4, find(1-tetSet)))'; zeros(4, 2 * newTets)']';

    index = numTets - newTets + 1;
    for i=1:numEdgesPicked
        if pickedEdges(i) == 0
            break
        end

        edge = pickedEdges(i);
        newEdge = newToOld(edge);
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
                newElements(:, index) = newToOld(elements(1:4, tet));
                newElements(endpoints(endpoint), index) = newEdge;
                index = index + 1;
            end            
        end
    end
    
    data = getTetDataRT(newElements', newNodes');
    newMaxVert = max(max(newElements));
    edgeToVert = containers.Map(cellstr(num2str(data.edges, "%d %d")), (newMaxVert+1):(newMaxVert+data.numEdges));
    edgeInd = [1, 2; 2, 3; 1, 3; 1, 4; 2, 4; 3, 4]';
    newEdges = sort(reshape(newElements(edgeInd, :), 2, []))';
    newEdges = values(edgeToVert, cellstr(num2str(newEdges, '%d %d')));
    newEdges = reshape(cell2mat(newEdges), 6, []);
    newElements = [newElements; newEdges];
    newNodes = [newNodes, reshape(mean(reshape(newNodes(:, data.edges), [], 2), 2), 3, [])];
    
    newModel = createpde('thermal', 'transient');
    try
        warning('off', 'pde:PDEModel:NodesNotInMeshWarnId');
        geometryFromMesh(newModel, newNodes, newElements);
        warning('on', 'pde:PDEModel:NodesNotInMeshWarnId');
    catch
        disp('hola');
    end
%     newQ = (buildIHh(newModel.Mesh)*newQ')';
    
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