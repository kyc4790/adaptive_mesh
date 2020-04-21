% model = createpde('thermal', 'transient');
% importGeometry(model, 'tetrahedron.stl');
% generateMesh(model, 'Hmin',15, 'GeometricOrder', 'quadratic');
% 
% mesh = buildMesh(model.Mesh);
% qOcta = MBO(mesh, OctaMBO, [], 1, 0);

function [lambdaH]=getLambdaH(mesh, qOcta, dg)
    Bflattened = getDerivOfF(mesh, qOcta);
    if nargin < 3 || isempty(dg)
        dg = getDerivOfG(mesh, qOcta);
    end
    lambdaH = Bflattened \ dg;
end
