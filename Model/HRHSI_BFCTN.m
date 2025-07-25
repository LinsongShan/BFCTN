function HRHSI = HRHSI_BFCTN(Org, MSI, HSI, P1, P2, P3, para)
% HRHSI_BFCTN - Hyperspectral and multispectral fusion via Bayesian Fully-Connected Tensor Network
%
% Inputs:
%   Org  - Ground truth tensor (high-res hyperspectral image)
%   MSI  - Observed high-resolution multispectral image
%   HSI  - Observed low-resolution hyperspectral image
%   P1,P2,P3 - Projection matrices (spatial and spectral degradation)
%   para - Structure containing algorithm parameters
%
% Output:
%   HRHSI - Fused high-resolution hyperspectral image

%% Size settings
[M, N, l] = size(MSI);      % MSI spatial & spectral size
[m, n, L] = size(HSI);      % HSI size
msi_patchsize = para.msi_patchsize; 
hsi_patchsize = para.msi_patchsize / para.ratio;

mparams.block_sz = [msi_patchsize, msi_patchsize]; 
hparams.block_sz = [hsi_patchsize, hsi_patchsize];
full.block_sz    = [msi_patchsize, msi_patchsize];

overlap    = para.overlap;
overlap_sz = [overlap, overlap];

mparams.sz = [M, N, l];
hparams.sz = [m, n, L];
full.sz    = [M, N, L];

num1 = (M - overlap) / (msi_patchsize - overlap);  % vertical blocks
num2 = (N - overlap) / (msi_patchsize - overlap);  % horizontal blocks

mparams.block_num = [num1, num2]; 
hparams.block_num = [num1, num2];
full.block_num    = [num1, num2];

%% Extract overlapping patches
MSI_blocks = My_ExtractBlocks(MSI, mparams, overlap_sz);                % Extract MSI patches
HSI_blocks = My_ExtractBlocks(HSI, hparams, overlap_sz / para.ratio);   % Extract HSI patches
O_blocks   = My_ExtractBlocks(Org, full, overlap_sz);                   % Extract GT patches

cluster = 1;

%% Process block by block along horizontal dimension
while cluster <= num2
    fprintf('*** Processing group %d of %d ***\n', cluster, num2);

    M_k = zeros(msi_patchsize, msi_patchsize, l, num1);    % MSI patches for this group
    H_k = zeros(hsi_patchsize, hsi_patchsize, L, num1);    % HSI patches
    O_k = zeros(msi_patchsize, msi_patchsize, L, num1);    % Ground truth patches

    for i = 1:num1
        idx = num1 * (cluster - 1) + i;
        M_k(:,:,:,i) = MSI_blocks(:,:,:,idx);
        H_k(:,:,:,i) = HSI_blocks(:,:,:,idx);
        O_k(:,:,:,i) = O_blocks(:,:,:,idx);
    end

    % Scaling to [0,255] range
    c = 255;
    M_k = M_k * c;
    H_k = H_k * c;
    O_k = O_k * c;

    % Hyperparameters for Bayesian inference
    MaxRank = para.MaxRank;
    MAX_Iter  = para.MAX_Iter;

    hyperparameters.lambda_a_0 = 1e-6;
    hyperparameters.lambda_b_0 = 1e-6;
    hyperparameters.alpha_c_0 = 1e-6;
    hyperparameters.alpha_d_0 = 1e-6;
    hyperparameters.beta_e_0  = 1e-6;
    hyperparameters.beta_f_0  = 1e-6;

    % Run BFCTN model
    VB = BFCTN_model(M_k, H_k, O_k, P1, P2, P3, MaxRank, hyperparameters);
    VB = VB.initialize();
    VB = VB.run(MAX_Iter);

    % Reconstruct fused result for this group
    ZZ = double(tnprod_new(VB.Factors));

    % Clip to valid range
    % max_v = max(O_k(:));
    % min_v = min(O_k(:));
    % ZZ(ZZ > max_v) = max_v;
    % ZZ(ZZ < min_v) = min_v;

    for i = 1:num1
        idx = num1 * (cluster - 1) + i;
        Z(:,:,:,idx) = ZZ(:,:,:,i);
    end 

    cluster = cluster + 1;
    fprintf('\n');
end

%% Merge blocks into full image and rescale
HRHSI = merge(Z, full, overlap_sz) / c;

end
