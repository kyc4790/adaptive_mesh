function [df]=getDerivOfF(mesh, qOcta)
    df=reshape((mesh.fineL*qOcta')', [], 1);
end