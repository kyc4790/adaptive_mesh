function overlay(info, first, second)
    length = size(info, 2);
    if(nargin < 3)
        first = 1;
        second = length - 1;
    end
        
    VisualizeResult(info(second).mesh, info(second).q);
    compareMesh(info(first).model, info(second).model, 0);
end