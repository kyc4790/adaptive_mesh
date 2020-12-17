function [qProj] = projBoundary(alpha, normals, qBdry)
% Projects qBdry onto SO3 with relaxed constraints 
%                   s.t. abs(W_n x - b) <= alpha
% 
% @param alpha      the threshold for relaxation
% @param normals    the normals of the boundary
% @param qBdry      the unprojected frame for the boundary

    numBdry = size(normals, 1);
    
    mat = OctaAlignMat(normals);
    W = multiprod(multitransp(mat), sparse(2:8, 1:7, ones(1, 7), 9, 7));
    
    % Define each of the 7 constraints to avoid tensors
    Wrows = repelem(1:numBdry, 9);
    Wcols = 1:(9*numBdry);
    const1 = sparse(Wrows, Wcols, reshape(W(:, 1, :), 1, []), numBdry, numBdry * 9);
    const2 = sparse(Wrows, Wcols, reshape(W(:, 2, :), 1, []), numBdry, numBdry * 9);
    const3 = sparse(Wrows, Wcols, reshape(W(:, 3, :), 1, []), numBdry, numBdry * 9);
    const4 = sparse(Wrows, Wcols, reshape(W(:, 4, :), 1, []), numBdry, numBdry * 9);
    const5 = sparse(Wrows, Wcols, reshape(W(:, 5, :), 1, []), numBdry, numBdry * 9);
    const6 = sparse(Wrows, Wcols, reshape(W(:, 6, :), 1, []), numBdry, numBdry * 9);
    const7 = sparse(Wrows, Wcols, reshape(W(:, 7, :), 1, []), numBdry, numBdry * 9);
    
    qBdryReshaped = reshape(qBdry, [], 1);
    sum_mag = sparse(repelem(1:numBdry, 9), 1:(9*numBdry), ones(1, 9*numBdry), numBdry, 9 * numBdry);
    
    cvx_begin quiet
        cvx_solver mosek
        cvx_precision high
        variable x(9*numBdry)
        minimize( norm( x - qBdryReshaped, 2) )
        subject to
            abs(const1 * x) <= alpha
            abs(const2 * x) <= alpha
            abs(const3 * x) <= alpha
            abs(const4 * x - sqrt(7/12)) <= alpha
            abs(const5 * x) <= alpha
            abs(const6 * x) <= alpha
            abs(const7 * x) <= alpha
            sum_mag * power(x, 2) <= 1
    cvx_end
    qProj = reshape(x, 9, []);
end