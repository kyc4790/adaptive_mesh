% plot(cell2mat({info.inner_iters}));
% title('inner iters');
meshes = {info.mesh};
[matSize1, matSize2, nv] = cellfun(@getSizeOfMatrixSolve, meshes);


yyaxis left
plot(1:5, cell2mat({info.inner_time}), (1:4) + 0.5, cell2mat({info.costTime}));
title('inner time vs cost time');
yyaxis right
plot(1:5, nv);


function [dim1, dim2, nv]=getSizeOfMatrixSolve(mesh)
    dim1 = 9*mesh.nv + size(mesh.bdryIdxFine, 1) * 8 + size( mesh.intIdxFine, 1) * 15;
    dim2 = dim1 * 2 - 9 * mesh.nv;
    nv = mesh.nv;
end

