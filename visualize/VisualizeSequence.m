function VisualizeSequence(info)
% Visualizes a sequence of edge splits. 
%
% @param info the output of adaptiveRemesh.
    axes = [];
    fig = figure;
    n = size(info, 2);

    for i=1:n
        meshData = info(i).mesh;
        q = info(i).q;
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
        ax1 = subplot(3, n, i);
        title(ax1, num2str(info(i).cost));
        try
            PlotIntegralCurves(frames, singularTetra, 'NumSeeds', 10000, 'Prune', true, 'ColorField', energy);
        catch ME
            disp(ME)
        end
        hold on;
        trisurf(bdry, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.05);
        view(3);
        axis image vis3d off;

        % Integral Curves
        ax2 = subplot(3, n, i+n);
        if(i < n)
            compareMesh(info(i).model, info(i+1).model, 0); 
        else
            compareMesh(info(1).model, info(n).model, 3);
        end
        tets = info(i).model.Mesh.Elements(1:4, :);
        vertices = info(i).model.Mesh.Nodes;
        T = triangulation(tets', vertices');
        [~, r] = incenter(T);
        edge_1 = vecnorm(vertices(:, (tets(1, :))) - vertices(:, (tets(2, :))));
        edge_2 = vecnorm(vertices(:, (tets(1, :))) - vertices(:, (tets(3, :))));
        edge_3 = vecnorm(vertices(:, (tets(1, :))) - vertices(:, (tets(4, :))));
        edge_4 = vecnorm(vertices(:, (tets(2, :))) - vertices(:, (tets(3, :))));
        edge_5 = vecnorm(vertices(:, (tets(2, :))) - vertices(:, (tets(4, :))));
        edge_6 = vecnorm(vertices(:, (tets(3, :))) - vertices(:, (tets(4, :))));
        tet_length = max([edge_1; edge_2; edge_3; edge_4; edge_5; edge_6]);
        tet_ratios = r./tet_length';
        [sorted, idx] = sort(tet_ratios);
        numSorted = size(sorted, 1);
        numSelected = int32(numSorted * 0.1);
        ax3 = subplot(3, n, i + 2*n);
        tetramesh(tets(:, idx(1:numSelected))', vertices' - mean(vertices'), tet_ratios(idx(1:numSelected), :), 'FaceAlpha',0.5, 'EdgeColor', 'none');
%         trisurf(bdry, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.05);
%         PlotIntegralCurves(frames, tetra, 'NumSeeds', 10000, 'Prune', true);
        view(3);
        axes = [axes, ax1, ax2, ax3];
        axis image vis3d off;

    end
    

    Link = linkprop(axes, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(fig, 'TheLink', Link);

    set(fig, 'color', 'white');
end

