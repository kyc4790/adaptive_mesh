plotData.bdry = mesh.bdry;
slider_plot(plotData);

function [] = slider_plot(plotData)
% Plot different plots according to slider location.
S.fh = figure('name','slider_plot');    
S.x = 0:.01:1;  % For plotting.         
S.LN = scatter3(plotData.x, plotData.y, plotData.z, 100, plotData.colors, 'filled');
axis equal
axis off
set(gcf, 'color', 'w');
S.plotData = plotData;
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 10 260 30],...
                 'min',0,'max', max(1.01*plotData.colors),'val',0,...
                 'sliderstep',[1/20 1/20],...
                 'callback',{@sl_call,S}); 
end
function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
plotData = S.plotData;
disp(h);
indexes = plotData.colors > h.Value;

scatter3(plotData.x(indexes, :), plotData.y(indexes, :), plotData.z(indexes, :), 100, plotData.colors(:, indexes), 'filled');
hold on
axis equal
axis off
set(gcf, 'color', 'w');
hold off
end