%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title  : Bayesian Fully-Connected Tensor Network for Hyperspectral-Multispectral Image Fusion
% Author : Linsong Shan, Laurence T. Yang, Zecan Yang, Changlong Li, Honglu Zhao, Xin Nie
% Date   : July 2025
%
% Description:
% This demo performs hyperspectral and multispectral image fusion using the proposed Bayesian Fully-Connected Tensor Network (BFCTN).
%
% Input:
%   - Simulated high-resolution multispectral image (HRMSI)
%   - Simulated low-resolution hyperspectral image (LRHSI)
%
% Output:
%   - Z: Fused high-resolution hyperspectral image
%   - Quantitative metrics: PSNR, RMSE, ERGAS, SAM, SSIM, UIQI, DD, CC
%   - Visualization of original and reconstructed images
%
% Dependencies:
%   - Model/: Model implementation
%   - Function/: Utility functions
%   - data/: Test data (e.g., chart_and_stuffed_toy_ms.mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

%% Add necessary paths
addpath(genpath('./Model'));
addpath(genpath('./data'));
addpath(genpath('./Function'));

%% Load data
load('chart_and_stuffed_toy_ms.mat');         % Load data
F = create_F();                                
orgTensor = msi(1:512, 1:512, :);              
orgTensor = (orgTensor - min(orgTensor(:))) / (max(orgTensor(:)) - min(orgTensor(:)));  % Normalize to [0,1]
szX = size(orgTensor);

% Visualize the original tensor in RGB
H_chanel = [28, 16, 1];                        % Band indices for visualization (e.g., RGB for CAVE)
orgImage = func_hyperImshow(orgTensor, H_chanel);
figure, imshow(orgImage, []), title('Original Tensor');

%% Define fusion scenario
sf = 4;                                        % Spatial downsampling factor
s0 = sf / 2;                                   
SNR = 35;                                      % Signal-to-noise ratio (dB)
degradation_type = 'average';  % Options: 'average', 'gaussian', 'motion'

%% Generate degraded observations
MSI = generate_HRMSI(orgTensor, F);
HSI = generate_LRHSI(orgTensor, sf, degradation_type);

% Add noise
MSI = Addnoise(MSI, SNR);
HSI = Addnoise(HSI, SNR);

%% Set model parameters
para.overlap       = 48;                      % Patch overlap size
para.msi_patchsize = 64;                      % MSI patch size

para.MaxRank = [0, 35,  4, 12;
                0,  0,  4, 12;
                0,  0,  0,  4;
                0,  0,  0,  0];               % FCTN ranks 

para.MAX_Iter  = 6;                           % Max VB iterations
para.ratio     = sf;                          

% Construct degradation kernel & projection
switch degradation_type
    case 'average'
        psf = fspecial('average', sf);
        P1 = creat_P(para.msi_patchsize, sf);
    case 'gaussian'
        psf_size = 4; sigma = 2;
        psf = fspecial('gaussian', psf_size, sigma);
        P1 = create_P_gaussian(para.msi_patchsize, sf, psf, psf_size);
    case 'motion'
        len = 4; theta = 30;
        psf = fspecial('motion', len, theta);
        P1 = create_P_motion(para.msi_patchsize, sf, psf);
    otherwise
        error('Unsupported degradation type');
end

otf = psf2otf(psf, szX(1:2));

%% Run BFCTN fusion
tic;
Z = HRHSI_BFCTN(orgTensor, MSI, HSI, P1, P1, F, para);
t = toc;

% Evaluate quality metrics
[psnr, rmse, ergas, sam, uiqi, ssim, DD, CC] = ...
    quality_assessment(double(im2uint8(orgTensor)), double(im2uint8(Z)), 0, 1.0/sf);

% Visualize fused result
ZImage = func_hyperImshow(Z, H_chanel);
figure, imshow(ZImage, []), title('Reconstructed Z');


