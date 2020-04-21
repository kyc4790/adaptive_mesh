function [d2f]=getHessianOfF(mesh, qOcta)
    d2f = sparse(9*mesh.nv, 9*mesh.nv);

    inds = ((1:mesh.nv)-1)*9;
    for i=1:9
         d2f(inds+i, inds+i) = mesh.fineL;
    end
end