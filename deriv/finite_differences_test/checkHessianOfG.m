function checkHessianOfG(mesh, qOcta)
    d2g = getHessianOfG(mesh);
    delta = 1e-6;
    epsilon = 1e-6;
    for j=1:size(qOcta, 2)
        for i=1:9
            newQPos = qOcta;
            newQNeg = qOcta;
            newQPos(i, j) = newQPos(i, j) + delta;
            newQNeg(i, j) = newQNeg(i, j) - delta;
            d2gPos = getDerivOfG(mesh, newQPos);
            d2gNeg = getDerivOfG(mesh, newQNeg);
            ind = (j-1)*9 + i;
            d2g_approx = (d2gPos - d2gNeg) / (2 * delta);
            real_d2g = d2g(((1:9*mesh.nv) + (ind - 1) * 9 * mesh.nv), :);
            diff = max(max(abs(real_d2g - d2g_approx)));
            if(diff > epsilon)
                disp([diff, i, j]);
            end
        end
    end
end