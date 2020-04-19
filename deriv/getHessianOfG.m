function [d2g]=getHessianOfG(mesh)
    mats = struct2cell(load('OctaMat.mat'));
    
    numBoundary = size(mesh.bdryIdx, 1);
    numInt = size(mesh.intIdx, 1);
    d2g = spalloc(9*mesh.nv*9*mesh.nv, numBoundary*8 + numInt*15, 9*numBoundary + 81*numInt*15);

    % first numBoundary constraints: |X|^2 - 1 => 2
    [rowIdx, toSet] = reformat(mesh.nv, mesh.bdryIdx, 2*eye(9));
    colIdx = reshape(repmat(1:numBoundary, 81, 1), [], 1);
    d2g(sub2ind(size(d2g), rowIdx, colIdx)) = toSet;
    % next constraints: W_n x - c => no second deriv

    % next contraints: x' P_j x => 2*C
    startInd = 8*numBoundary;
    numInt = size(mesh.intIdx, 1);

    for i=1:15
        P = mats{i};
        C = P(2:10, 2:10);
        [rowIdx, toSet] = reformat(mesh.nv, mesh.intIdx, C);
        colIdx = startInd + reshape(repmat(1:numInt, 81, 1), [], 1);
        d2g(sub2ind(size(d2g), rowIdx, colIdx)) = toSet;
        startInd = startInd + numInt;
    end

    function [idx, toSet]=reformat(nv, inds, mat)
        % say inds is 2. then you need 10-18, nv*
        num = repmat(1:9, 9, 1);
        preTranspose = nv*9*(num' - 1) + num;
        numInds = size(inds, 1);
        preTranspose = repmat(preTranspose', 1, numInds);
        indsRep = repmat(inds, 1, 81);
        indsRep = (reshape(indsRep', 9, [])-1)*(nv*9+1);
        idx = reshape(indsRep + preTranspose, [], 1);
        
        toSet=repmat(reshape(mat', 1, []), 1, numInds);
    end
end
