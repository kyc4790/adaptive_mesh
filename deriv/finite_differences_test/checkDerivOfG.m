function checkDerivOfG(mesh, qOcta)
    dg = getDerivOfG(mesh, qOcta);
    delta = 1e-6;
    epsilon = 1e-6;
    for i=1:9
        for j=1:size(qOcta, 2)
            newQPos = qOcta;
            newQNeg = qOcta;
            newQPos(i, j) = newQPos(i, j) + delta;
            newQNeg(i, j) = newQNeg(i, j) - delta;
            gPos = getG(mesh, newQPos);
            gNeg = getG(mesh, newQNeg);
            g_diff = gPos - gNeg;
            dg_approx = g_diff / (2*delta);
            dg_real = dg(j*9 - 9 + i, :)';
            diff = max(max(abs(dg_approx - dg_real)));
            if(diff > epsilon)
                disp([i, j, diff]);
            end
        end
    end
end