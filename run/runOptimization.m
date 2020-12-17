% Runs the optimization for a specified mesh.

res = [3];
iters = [5];
filename = 'notch5_vtk';
for i=1:size(res, 2)
    tic;
%     model = createpde('thermal', 'transient');
%     importGeometry(model, strcat(filename, '.stl'));
%     generateMesh(model, 'Hmax', res);

    model = assignModel(getModel(LoadNotch5()));

    [qOcta, newModel, info] = adaptiveRemesh(model, [], true, 1e-3, iters(i));
%     newMesh = buildMesh(newModel.Mesh);
%     VisualizeResult(newMesh, qOcta);
%     compareMesh(model, newModel);
    save(strcat('octa_', filename, '_', int2str(res(i)), '.mat'), 'info');
    % visualizeMesh(newModel.Mesh.Nodes, newModel.Mesh.Elements);
end
