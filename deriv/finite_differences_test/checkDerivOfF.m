function checkDerivOfF(mesh, qOcta)
    df = getDerivOfF(mesh, qOcta);
    delta = 1e-6;
    epsilon = 1e-6;
    for i=1:9
        for j=1:size(qOcta, 2)
            newQPos = qOcta;
            newQNeg = qOcta;
            newQPos(i, j) = newQPos(i, j) + delta;
            newQNeg(i, j) = newQNeg(i, j) - delta;
            fPos = getF(mesh, newQPos);
            fNeg = getF(mesh, newQNeg);
            f_diff = fPos - fNeg;
            df_approx = f_diff / (2*delta);
            df_real = df(j*9 - 9 + i, :)';
            diff = max(max(abs(df_approx - df_real)));
            if(diff > epsilon)
                disp([i, j, diff]);
            end
        end
    end
end