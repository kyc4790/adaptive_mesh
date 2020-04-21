function [f]=getF(mesh, qOcta)
    f = 0.5 * sum(dot(qOcta' , mesh.fineL * qOcta'));
end