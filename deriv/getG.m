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
    mat = OctaAlignMat(mesh.bdryNormals);
    W = multiprod(multitransp(mat), sparse(2:8, 1:7, ones(1, 7), 9, 7));
    res = multiprod(multitransp(W), q(1:9, mesh.bdryIdx));
    q_inds = repelem(1:numBdry, 7, 1);
    const_inds = repelem(1:7, numBdry, 1)'; 
    inds = sub2ind(size(res), const_inds, q_inds, q_inds);
    b = [0 0 0 sqrt(7/12) 0 0 0]';
    B = repmat(b, 1, numBdry);
    vals = res(inds) - B;
    
    for i=1:7
        g(numBdry + i + (1:numBdry)*7 - 7) = vals(i, :);
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