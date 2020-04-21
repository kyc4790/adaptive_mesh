model = createpde('thermal', 'transient');
importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmin',35, 'GeometricOrder', 'quadratic');

qOcta = [];
old_cost = -1;
epsilon = 1e-2;

tic;
for i=1:20
    mesh = buildMesh(model.Mesh);

    [qOcta, ~, info] = MBO(mesh, OctaMBO, qOcta, 1, 0);
    cost = info(length(info)).cost; 
    cost = cost / length(qOcta);
    if abs((old_cost - cost)/ old_cost) < epsilon
        break;
    end
    [~, costs, qOctaFine] = getLambda(model, mesh, qOcta);
    [model, qOcta, newNodes, newElements] = edgeSplitXT(model, costs, qOctaFine);
    old_cost = cost;
    fprintf("t = %3.3gs, average_cost = %3.6g, inner_iters=%d\n", toc, cost, length(info));
end

