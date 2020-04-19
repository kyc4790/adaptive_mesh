% model = createpde('thermal', 'transient');
% importGeometry(model, 'tetrahedron.stl');
% generateMesh(model, 'Hmin',15, 'GeometricOrder', 'quadratic');
% 
% mesh = buildMesh(model.Mesh);
% qOcta = MBO(mesh, OctaMBO, [], 1, 0);

function [lambdaH]=getLambdaH(mesh, qOcta, dg)
    Bflattened = 2*reshape((mesh.fineL*qOcta')', [], 1);
    if nargin < 3 || isempty(dg)
        dg = getDerivOfG(mesh, qOcta);
    end
    lambdaH = Bflattened \ dg;
end
