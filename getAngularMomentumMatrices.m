function [Lx, Ly, Lz]=getAngularMomentumMatrices()
    L = zeros(9, 9);
    L(8, 1) = sqrt(2);
    L(7, 2) = sqrt(7/2);
    L(6, 3) = 3/sqrt(2);

    Lx = fliplr(fliplr(L)') + L;
    Lx(5, 4) = sqrt(10);
    Lx = Lx - Lx';

    Ly = -flipud(L)+fliplr(L)';
    Ly(6, 5) = sqrt(10);
    Ly = Ly - Ly';

    Lz = zeros(9, 9);
    Lz(1, 9) = 4;
    Lz(2, 8) = 3;
    Lz(3, 7) = 2;
    Lz(4, 5) = 1;
    Lz = Lz - Lz';
end