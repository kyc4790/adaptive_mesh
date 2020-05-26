function compareMesh(model, newModel, numTimesSplit)
    if(nargin < 3)
        numTimesSplit = 1;
    end
    
    numOldNodes = max(max(model.Mesh.Elements(1:4, :)));
    numNewNodes = max(max(newModel.Mesh.Elements(1:4, :)));

    elements = newModel.Mesh.Elements(1:4, :);
    numTets = size(elements, 2);

    elem = reshape(elements, [], 1);
    elem = repmat(elem, 1, numNewNodes - numOldNodes);
    inds = repmat(numOldNodes + 1:numNewNodes, numTets * 4, 1);
    desired = (elem == inds);

    desired = sum(double(desired), 2);
    desired = reshape(desired, 4, []);
    desired = sum(desired);

    visualizeMesh(newModel.Mesh.Nodes(:, 1:numNewNodes), elements(:, desired > numTimesSplit));
end
