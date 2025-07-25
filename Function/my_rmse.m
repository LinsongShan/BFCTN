function [rmse] = my_rmse(ground_truth, estimated, ignore_edges)

% Ignore borders
y = ground_truth(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);
x = estimated(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);

% Size, bands, samples 
sz_x = size(x);
n_bands = sz_x(3);
n_samples = sz_x(1)*sz_x(2);

% RMSE
aux = sum(sum((x - y).^2, 1), 2)/n_samples;
rmse = sqrt(sum(aux, 3)/n_bands);


