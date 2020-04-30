function [d2f]=getHessianOfF(mesh, qOcta)
    d2f = sparse(9*mesh.newnv, 9*mesh.newnv);

    inds = ((1:mesh.newnv)-1)*9;
    for i=1:9
         d2f(inds+i, inds+i) = mesh.fineL;
    end
end