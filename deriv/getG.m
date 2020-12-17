function [g]=getG(mesh, q)
    inds = 7:15;
    if(size(q, 1) == 9)
        inds = 1:9;
    end
    mats = struct2cell(load('OdecoMatSph.mat'));
    qy = q;
    if(size(q, 1) == 9)
       qy = [ones(1, size(q, 2)); q];
       mats = LoadOctaMatsScaled; % struct2cell(load('OctaMat.mat'));
    end
    numBdry = size(mesh.bdryIdx, 1);
    numInt = size(mesh.intIdx, 1);
    g = zeros(numBdry * 8 + numInt * size(mats, 1), 1);
    
    % |x|^2 - 1
    g(1:numBdry) = diag(q(inds, mesh.bdryIdx)' * q(inds, mesh.bdryIdx)) - 1;
    
    %W_n x - constant
    normalVecs = num2cell(mesh.bdryNormals, 2);
    W = cellfun(@(x) getDot(x), normalVecs, 'UniformOutput', false);
    wnx = cell2mat(W)*q(inds, mesh.bdryIdx);
    constant = zeros(7, 1);
    constant(4) = sqrt(7/12);
    for i=1:7
        g(numBdry + i + (1:numBdry)*7 - 7) = wnx(sub2ind(size(wnx), 7*(1:numBdry) - 7 + i, 1:numBdry)) - constant(i);
    end
    
    % x^TPx


    qy = qy(:, mesh.intIdx);
    startInd =  numBdry * 8;
    for i=1:size(mats, 1)
        P = mats{i};
        mat = qy'*P*qy;
        g(startInd + 1:startInd + numInt) = diag(mat);
        startInd = startInd + numInt;
    end
    
    function [result]=getDot(normalVec)
        [Lx, Ly, Lz] = getAngularMomentumMatrices();
        result = expm(-1*(normalVec(1) * Lx + normalVec(2) * Ly + normalVec(3) * Lz));

        result = result(2:8, :);
    end
    
end