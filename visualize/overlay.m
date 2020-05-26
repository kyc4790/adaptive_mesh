function overlay(info)
    length = size(info, 2);
    VisualizeResult(info(length).mesh, info(length).q);
    compareMesh(info(1).model, info(length-1).model);
end