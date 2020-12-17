function visualizeVariation(weightedMesh, info)
%VISUALIZEVARIATION Plays a movie of singular curves over time
%   @param info             the output of runVariation
%   @param weightedMesh     the MBOmesh associated with the tet mesh.
    n = size(info, 2);
    pause on;
    while(true)
        for index=1:n
            for j=1:size(info(index).i, 2)
                render(weightedMesh, info(index).i(j).q);
                pause(0.5);
            end
        end
    end

%     fig = uifigure();
%     
%     display(n);
% 
%     sld = uislider(fig,...
%                'Position',[100 75 120 3],...
%                'Limits', [1, n],...
%                'Value', 1,...
%                'ValueChangingFcn',@(sld,event) sliderMoving(sld,mesh, weightedMesh, info1, info2));
end

function sliderMoving(sld, mesh, weightedMesh, info1, info2)
    index = sld.Value;
    if(index <= size(info1, 2))
        render(mesh, info1(index).q);
    else
        render(weightedMesh, info2(index - size(info1, 2)).q);
    end
end

function render(meshData, q)
    cla;
    frames = Coeff2Frames(q);
    tetra = meshData.tetra;
    bdry = meshData.bdry;

    if size(q, 1) == 9 % Octahedral
        energy = dot(q.', meshData.L * q.', 2);
    else % Odeco
        energy = dot(q(7:15, :).', meshData.L * q(7:15, :).', 2);
    end


    frames = frames ./ vecnorm(frames, 2, 1);

    % Singularities only
    % singularIdx = find(energy > mean(energy) + std(energy));
    % singularTets = vertexAttachments(tetra, singularIdx)';
    % singularTets = unique(cell2mat(singularTets)');
    singularTets = ExtractSingularities(frames, tetra);
    warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    singularTetra = triangulation(tetra(singularTets, :), tetra.Points);
    warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');
    PlotIntegralCurves(frames, singularTetra, 'NumSeeds', 10000, 'Prune', true, 'ColorField', energy);
    hold on;
    trisurf(bdry, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.05);
    view(3);
    axis image vis3d off;
end

