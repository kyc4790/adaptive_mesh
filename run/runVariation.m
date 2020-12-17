% Runs the original altering of weights by changing the function
% info => contains all the information throughout the run
% model => the mesh (triangulation) used
% weightedMesh => the mesh with additional details (including weights)

filename = 'penguin';

model = getModel(loadMesh(filename));
weightedMesh = buildWeightedMesh(model);
weightMax = max(weightedMesh.weights);

info = struct;
q = [];
numIters = 4;
numSecondIters = 2;
for i=1:numIters
    tic;
    [q, ~, info1] = MBO(adjustWeights(model, weightedMesh, @(x) x.^(i-1)), OctaMBO, q, 1, 0, true, 1000);
    info(i).i = info1;
    fprintf('%d %f\n', i, toc);
end

m = max(max(weightedMesh.weights .^(numIters - 1)));
b = numIters * max((weightedMesh.weights + 1).*log(1 + weightedMesh.weights));

for i=numIters + 1:(numIters + numSecondIters)
    tic;
    [q, ~, info1] = MBO(adjustWeights(model, weightedMesh, @(x) (log(1+x).^(i - numIters + b))), OctaMBO, q, 1, 0, true, 2);
    info(i).i = info1;
    fprintf('%d %f\n', i, toc);
end

save(strcat(filename, '_variation_poly_log_power.mat'), 'info', 'model', 'weightedMesh');

function [mesh] = loadMesh(filename)
    if(strcmp(filename, 'penguin'))
        mesh = LoadPenguinMesh();
    elseif(strcmp(filename, 'die'))
        mesh = LoadDieMesh();
    elseif(strcmp(filename, 'candle'))
        mesh = LoadCandleMesh();
    end
end