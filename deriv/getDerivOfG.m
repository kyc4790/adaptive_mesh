function [dg]=getDerivOfG(mesh, qOcta)
    mats = struct2cell(load('OctaMat.mat'));
    
    qFlattened = reshape(qOcta, 1, []);
    numBoundary = size(mesh.bdryIdx, 1);
    numInt = size(mesh.intIdx, 1);
    nv = mesh.nv;
    dg = spalloc(9*nv, numBoundary*8 + numInt*15, 9*(numBoundary*8 + numInt*15));

    % first numBoundary constraints: |X|^2 - 1
    rowIdx = reshape((repmat((mesh.bdryIdx-1)*9, 1, 9) + repmat(1:9, numBoundary, 1))', [], 1);
    colIdx = reshape(repmat(1:numBoundary, 9, 1), [], 1);
    dg(sub2ind(size(dg), rowIdx, colIdx)) = 2*qFlattened(rowIdx);

    % next constraints: W_n x - c
    normalVecs = num2cell(mesh.bdryNormals, 2);
    W = cellfun(@(x) getDot(x), normalVecs, 'UniformOutput', false);
    WMat = cell2mat(W);
    WRowIdx = reshape(repmat(rowIdx, 1, 7)', [], 1)';
    WColIdx = numBoundary + reshape(repmat(reshape(repmat(((1:numBoundary) - 1) * 7, 9, 1), [], 1), 1, 7)' + repmat(repmat(1:7, 9, 1)', 1, numBoundary), 1, []);
    WMatInd = repmat(1:63, 1, numBoundary)';
    WInd = reshape(repmat(1:numBoundary, 63, 1), [], 1);
    dg(sub2ind(size(dg), WRowIdx, WColIdx)) = WMat(sub2ind(size(WMat), WInd, WMatInd));

    % next contraints: x' P_j x
    startInd = 8*numBoundary;

    for i=1:15
        P = mats{i};
        B = P(2:10, 1);
        C = P(2:10, 2:10);
        rowIdx = reshape((repmat((mesh.intIdx-1)*9, 1, 9) + repmat(1:9, numInt, 1))', [], 1);
        colIdx = startInd + reshape(repmat(1:numInt, 9, 1), [], 1);
        cellQ = num2cell(qOcta(:,mesh.intIdx), 1);
        deriv = cellfun(@(x) getDeriv(B, C, x), cellQ, 'UniformOutput', false);
        derivFlattened = reshape(cell2mat(deriv), [], 1);
        dg(sub2ind(size(dg), rowIdx, colIdx)) = derivFlattened;
        startInd = startInd + numInt;
    end

    function [result]=getDeriv(B, C, x) 
        result = 2*B + 2*C*x;
    end

    function [result]=getDot(normalVec)
        [Lx, Ly, Lz] = getAngularMomentumMatrices();
        result = expm(normalVec(1) * Lx + normalVec(2) * Ly + normalVec(3) * Lz);
        result = result(2:8, :);
        result = reshape(result, 1, []);
    end
end
