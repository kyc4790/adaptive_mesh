clear; close all; clc;

model = createpde('thermal', 'transient');
importGeometry(model, 'tetrahedron.stl');
generateMesh(model, 'Hmin', 1, 'GeometricOrder', 'quadratic');

IHh = buildIHh(model.Mesh);
nodes = model.Mesh.Nodes;
elements = model.Mesh.Elements;
nodes = nodes - mean(nodes, 2);

temp = size(nodes);
maxNodes = temp(1, 2);
maxVertex = max(elements(1:4, :), [], 'all');

mesh = buildMesh(model.Mesh);


qOcta = MBO(mesh, OctaMBO, [], 1, 0);
qOcta = OctaManopt(mesh, qOcta);

A = [eye(maxNodes); IHh'];
B = [2*mesh.fineL*IHh*qOcta'; zeros(maxVertex, 9)];

lambda = B\A;
norms = vecnorm(lambda);

x = nodes(1, :)';
y = nodes(2, :)';
z = nodes(3, :)';
plotData = struct();
plotData.x = x;
plotData.y = y;
plotData.z = z;
plotData.colors = norms + eps;
scatter3(x, y, z, 100, norms + eps, 'filled');
slider_plot(plotData);
hold on
axis equal
axis off
set(gcf, 'color', 'w');
hold off

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
S.LN = plot(S.x,S.x,'r');        
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 10 260 30],...
                 'min',1,'max',10,'val',1,...
                 'sliderstep',[1/20 1/20],...
                 'callback',{@sl_call,S}); 
end
function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
set(S.LN,'ydata',S.x.^get(h,'value'));
end