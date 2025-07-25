function [psnr,rmse, ergas, sam, uiqi,ssim,DD,CCS] = quality_assessment(ground_truth, estimated, ignore_edges, ratio_ergas)

% Ignore borders
y = ground_truth(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);
x = estimated(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);

% Size, bands, samples 
sz_x = size(x);
n_bands = sz_x(3);
n_samples = sz_x(1)*sz_x(2);

% RMSE
aux = sum(sum((x - y).^2, 1), 2)/n_samples;
rmse_per_band = sqrt(aux);
rmse = sqrt(sum(aux, 3)/n_bands);

% ERGAS
mean_y = sum(sum(y, 1), 2)/n_samples;
ergas = 100*ratio_ergas*sqrt(sum((rmse_per_band ./ mean_y).^2)/n_bands);


q_band = zeros(1, n_bands);
for idx1=1:n_bands
    q_band(idx1)=img_qi(ground_truth(:,:,idx1), estimated(:,:,idx1), 32);
end

ssim=cal_ssim(ground_truth, estimated,0,0);
psnr=csnr(ground_truth, estimated,0,0);