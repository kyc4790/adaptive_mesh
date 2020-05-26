function [qOcta, model, info]=adaptiveRemesh(model, qOcta, saveIterates, epsilon, addInfo)
    if(nargin < 2)
        qOcta = [];
    end
    
    if(nargin < 3)
        saveIterates = false;
    end
    
    if(nargin < 4)
        epsilon = 1e-3;
    end
    
    start = 0;
    old_cost = -1;

    % Load from previous run
    if(nargin == 5)
        info = addInfo;
        start = length(info);
        model = info(start - 1).model;
        mesh = buildMesh(model.Mesh);
        
        tic;
        [costs, qOctaFine] = getDerivConstraint(model, mesh, qOcta);
        info(start).costTime = toc;
        tic;
        [model, qOcta] = edgeSplitXT(model, costs, qOctaFine);
        info(start).edgeSplitTime = toc;
        if(saveIterates)
            info(start).model = model;
            info(start).costs = costs;
        end
        fprintf("cost_t = %3.3gs, edge_t = %3.6gs\n", info(start).costTime, info(start).edgeSplitTime);
    end
    
    numIter = 18;

    for i=start+1:numIter
        tic;
        mesh = buildMesh(model.Mesh);

        [qOcta, ~, inner_info] = MBO(mesh, OctaMBO, qOcta, 1, 0);
        cost = inner_info(length(inner_info)).cost; 
%         cost = cost;
        
        info(i).inner_iters = length(inner_info);
        info(i).inner_time = toc;
        info(i).cost = cost;
        if(saveIterates)
            info(i).q = qOcta;
            info(i).mesh = mesh;
        end
        
        fprintf("t = %3.3gs, cost = %3.6g, inner_iters=%d\n", info(i).inner_time, cost, length(inner_info));

        
        if abs((old_cost - cost)/ old_cost) < epsilon
            fprintf("exited via break\n");
            break;
        end
        if(i < numIter)
            tic;
            [costs, qOctaFine] = getDerivConstraint(model, mesh, qOcta);
            info(i).costTime = toc;
            tic;
            [model, qOcta] = edgeSplitXT(model, costs, qOctaFine);
            info(i).edgeSplitTime = toc;
            if(saveIterates)
                info(i).model = model;
                info(i).costs = costs;
            end
            fprintf("cost_t = %3.3gs, edge_t = %3.6gs\n", info(i).costTime, info(i).edgeSplitTime);
            old_cost = cost;
        end
    end
end