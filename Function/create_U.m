function [U] = create_U(org_bands, ratio)
%% Generate spectral downsampling matrix
% Author: 
% Time: 2024-04-03
% Input:
%    org_band - The spectral bands of a hyperspectral image HRHSI
%    ratio - The downsampling ratio of the spectrum
% Output:
%    U - spectral downsampling matrix 
%% Main Function

    bands = ceil(org_bands/ratio);
    U = zeros(bands, org_bands);
    j = 1;
    for i = 1 : bands
        U(i, j) = 1;
        j = j + ratio;
    end

end