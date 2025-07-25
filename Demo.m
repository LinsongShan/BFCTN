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
load('chart_and_stuffed_toy_ms.mat');         % Load MSI data
F = create_F();                                % Create spectral response matrix
orgTensor = msi(1:512, 1:512, :);              % Crop to 512×512 spatial region
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

% Generate degraded observations
MSI = generate_HRMSI(orgTensor, F);           % High-resolution MSI
HSI = generate_LRHSI(orgTensor, sf);          % Low-resolution HSI
% HSI = generate_LRHSI_2(orgTensor, sf);      % Optional: Gaussian-based downsampling

% Add Gaussian noise to observations
MSI = Addnoise(MSI, SNR);
HSI = Addnoise(HSI, SNR);

% Define PSF and convert to OTF
psf = fspecial('average', sf);                % Box filter PSF
% psf = fspecial('gaussian', 4, 2);           % Optional: Gaussian PSF
otf = psf2otf(psf, szX(1:2));                 % Optical Transfer Function

%% Set model parameters
para.overlap       = 48;                      % Patch overlap size
para.msi_patchsize = 64;                      % MSI patch size

para.MaxRank = [0, 35,  4, 12;
                0,  0,  4, 12;
                0,  0,  0,  4;
                0,  0,  0,  0];               % Tensor rank structure

para.MAX_Iter  = 6;                           % Max VB iterations
para.ratio     = sf;                          % Upsampling ratio

% Create projection matrices
P1 = creat_P(para.msi_patchsize, sf);         % Spatial degradation matrix
P3 = F;                                        % Spectral response matrix

%% Run BFCTN fusion
tic;
Z = HRHSI_BFCTN(orgTensor, MSI, HSI, P1, P1, P3, para);
t = toc;

% Evaluate quality metrics
[psnr, rmse, ergas, sam, uiqi, ssim, DD, CC] = ...
    quality_assessment(double(im2uint8(orgTensor)), double(im2uint8(Z)), 0, 1.0/sf);

% Visualize fused result
ZImage = func_hyperImshow(Z, H_chanel);
figure, imshow(ZImage, []), title('Reconstructed Z');


