function checkHessianOfF(mesh, qOcta)
    d2f = getHessianOfF(mesh);
    delta = 1e-6;
    epsilon = 1e-6;
    for j=1:size(qOcta, 2)
        for i=1:9
            newQPos = qOcta;
            newQNeg = qOcta;
            newQPos(i, j) = newQPos(i, j) + delta;
            newQNeg(i, j) = newQNeg(i, j) - delta;
            d2fPos = getDerivOfF(mesh, newQPos);
            d2fNeg = getDerivOfF(mesh, newQNeg);
            ind = (j-1)*9 + i;
            d2f_approx = (d2fPos - d2fNeg) / (2 * delta);
            real_d2f = d2f(ind, :)';
            diff = max(max(abs(real_d2f - d2f_approx)));
            if(diff > epsilon)
                disp([diff, i, j]);
            end
        end
    end
end