plotData.bdry = mesh.bdry;
slider_plot(plotData);

function [] = slider_plot(plotData)
% Plot different plots according to slider location.
S.fh = figure('units','pixels',...
              'position',[300 300 300 300],...
              'menubar','none',...
              'name','slider_plot',...
              'numbertitle','off',...
              'resize','off');    
S.x = 0:.01:1;  % For plotting.         
S.ax = axes('unit','pix',...
            'position',[20 80 260 210]);
S.LN = scatter3(plotData.x, plotData.y, plotData.z, 100, plotData.colors, 'filled');
S.plotData = plotData;
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 10 260 30],...
                 'min',0,'max', 0.5,'val',0,...
                 'sliderstep',[1/20 1/20],...
                 'callback',{@sl_call,S}); 
end
function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
plotData = S.plotData;
disp(h);
indexes = plotData.colors > h.Value;

trimesh(plotData.bdry, 'FaceAlpha', 0, 'EdgeColor', 'k');
hold on
scatter3(plotData.x(indexes, :), plotData.y(indexes, :), plotData.z(indexes, :), 100, plotData.colors(:, indexes));
axis equal
axis off
set(gcf, 'color', 'w');
hold off
% text(plotData.x(indexes, :), plotData.y(indexes, :), plotData.z(indexes, :), num2str(plotData.colors(indexes, :)));
end