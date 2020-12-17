res = [2, 1];
filename = 'cylinder';
for i=1:size(res, 2)
    fprintf('%f\n', res(i));
    model = createpde('thermal', 'transient');
    importGeometry(model, strcat(filename,'.stl'));
    generateMesh(model, 'Hmax', res(i), 'GeometricOrder', 'quadratic');

    mesh = buildMesh(model.Mesh, false);
    [qOcta, ~, info] = MBO(mesh, OctaMBO, [], 1, 0, true);

    save(strcat('octa_no_split_', filename, '_', int2str(res(i)), '.mat'), 'info', 'mesh');
    % visualizeMesh(newModel.Mesh.Nodes, newModel.Mesh.Elements);
end