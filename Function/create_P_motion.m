% function S  = create_P_motion(patchsize, space_ratio, psf, psf_size)
%       [U, S, V] = svd(psf);
%         u1 = U(:, 1); % 高度方向近似向量
%         v1 = V(:, 1)'; % 宽度方向近似向量
% 
%         % 归一化
%         u1 = u1 / sum(u1);
%         v1 = v1 / sum(v1);
% 
%         % 构造降采样矩阵 S (高度方向)
%         H_original = patchsize;
%         H_lr = H_original / space_ratio;
%         S = zeros(H_lr, H_original);
%         for i = 1:H_lr
%             center = (i-1)*space_ratio + (space_ratio/2);
%             start = max(1, round(center - (psf_size-1)/2));
%             stop = min(H_original, round(center + (psf_size-1)/2));
%             S(i, start:stop) = u1(1:(stop-start+1));
%         end
% 
%         % % 构造降采样矩阵 T (宽度方向)
%         % W_original = 256;
%         % W_lr = W_original / sf;
%         % T = zeros(W_original, W_lr);
%         % for j = 1:W_lr
%         %     center = (j-1)*sf + (sf/2);
%         %     start = max(1, round(center - (len-1)/2));
%         %     stop = min(W_original, round(center + (len-1)/2));
%         %     T(start:stop, j) = v1(1:(stop-start+1))';
%         % end
% end


% function S = create_P_motion(patchsize, sf, psf)
%     % 输入参数:
%     %   patchsize : 原始图像尺寸（假设为正方形，如256）
%     %   sf        : 降采样比例（如4）
%     %   psf       : 运动模糊核（由 fspecial('motion', len, theta) 生成）
% 
%     % --- 步骤1：将运动模糊核转为方阵 ---
%     [h, w] = size(psf);
%     len = max(h, w); % 运动模糊核的理论长度（如4）
%     psf_square = zeros(len, len);
%     center = ceil(len/2);
% 
%     % 自适应填充（兼容任意psf尺寸）
%     if h == 1 % 水平运动模糊（1×len）
%         psf_square(center, :) = psf;
%     elseif w == 1 % 垂直运动模糊（len×1）
%         psf_square(:, center) = psf;
%     else % 非对称模糊核（如3×5），强制投影到中心行/列
%         row_proj = sum(psf, 1); % 行投影（1×w）
%         col_proj = sum(psf, 2); % 列投影（h×1）
%         if sum(row_proj) >= sum(col_proj)
%             psf_square(center, :) = row_proj / sum(row_proj);
%         else
%             psf_square(:, center) = col_proj / sum(col_proj);
%         end
%     end
% 
%     % --- 步骤2：SVD分解得到可分离向量 ---
%     [U, ~, V] = svd(psf_square);
%     u1 = U(:, 1); % 高度方向向量
%     v1 = V(:, 1)'; % 宽度方向向量
% 
%     % 归一化（保证能量守恒）
%     u1 = u1 / sum(u1);
%     v1 = v1 / sum(v1);
% 
%     % --- 步骤3：构造降采样矩阵 S（高度方向）---
%     H_original = patchsize;
%     H_lr = floor(H_original / sf);
%     S = zeros(H_lr, H_original);
% 
%     for i = 1:H_lr
%         center_pos = (i-1)*sf + ceil(sf/2);
%         start = max(1, center_pos - floor(len/2));
%         stop = min(H_original, center_pos + floor(len/2));
%         valid_range = start:stop;
%         u1_segment = u1(1:numel(valid_range));
%         S(i, valid_range) = u1_segment / sum(u1_segment);
%     end
% 
%     % --- 步骤4：构造降采样矩阵 T（宽度方向）---
%     W_original = patchsize;
%     W_lr = floor(W_original / sf);
%     T = zeros(W_original, W_lr);
% 
%     for j = 1:W_lr
%         center_pos = (j-1)*sf + ceil(sf/2);
%         start = max(1, center_pos - floor(len/2));
%         stop = min(W_original, center_pos + floor(len/2));
%         valid_range = start:stop;
%         v1_segment = v1(1:numel(valid_range))';
%         T(valid_range, j) = v1_segment / sum(v1_segment);
%     end
% end


function S = create_P_motion(patchsize, sf, psf)
    % 输入参数:
    %   patchsize : 原始图像尺寸（假设为正方形，如256）
    %   sf        : 降采样比例（如4）
    %   psf       : 运动模糊核（由 fspecial('motion', len, theta) 生成）
    
    % --- 步骤1：运动模糊核的对称化处理 ---
    [h, w] = size(psf);
    if h == 1 || w == 1 % 向量型运动模糊核
        % 构造对称核以保留方向信息
        kernel_length = max(h, w);
        symmetric_psf = zeros(kernel_length, kernel_length);
        center = ceil(kernel_length/2);
        if h > w % 垂直运动模糊
            symmetric_psf(:, center) = psf;
        else % 水平运动模糊
            symmetric_psf(center, :) = psf;
        end
    else % 非对称矩阵（如斜向运动模糊）
        % 双线性插值扩展为方阵
        [x, y] = meshgrid(1:w, 1:h);
        [xi, yi] = meshgrid(linspace(1, w, max(h, w)), linspace(1, h, max(h, w)));
        symmetric_psf = interp2(x, y, psf, xi, yi, 'linear', 0);
    end
    
    % --- 步骤2：改进的SVD分解策略 ---
    [U, S_vec, V] = svd(symmetric_psf);
    energy_ratio = S_vec(1,1) / sum(diag(S_vec)); % 第一奇异值能量占比
    
    % 选择主导方向（保留更多方向信息）
    if energy_ratio < 0.85 % 能量分散时保留前两主成分
        u1 = U(:, 1) * sqrt(S_vec(1,1)) + U(:, 2) * sqrt(S_vec(2,2));
        v1 = (V(:, 1) * sqrt(S_vec(1,1)) + V(:, 2) * sqrt(S_vec(2,2)))';
    else
        u1 = U(:, 1) * sqrt(S_vec(1,1));
        v1 = V(:, 1)' * sqrt(S_vec(1,1));
    end
    
    % 归一化处理
    u1 = u1 / sum(u1);
    v1 = v1 / sum(v1);
    
    % --- 步骤3：自适应降采样矩阵构造 ---
    H_original = patchsize;
    H_lr = floor(H_original / sf);
    S = zeros(H_lr, H_original);
    
    % 构造高度方向矩阵
    kernel_radius = floor(length(u1)/2);
    for i = 1:H_lr
        center_pos = round((i-0.5)*sf); % 改进的中心位置计算
        start = max(1, center_pos - kernel_radius);
        stop = min(H_original, center_pos + kernel_radius);
        
        % 动态截取核的对应部分
        kernel_start = max(1, kernel_radius + 2 - (center_pos - start));
        kernel_stop = min(length(u1), kernel_radius + 1 + (stop - center_pos));
        valid_kernel = u1(kernel_start:kernel_stop);
        
        S(i, start:stop) = valid_kernel / sum(valid_kernel);
    end
    
    % 构造宽度方向矩阵
    W_original = patchsize;
    W_lr = floor(W_original / sf);
    T = zeros(W_original, W_lr);
    
    kernel_radius = floor(length(v1)/2);
    for j = 1:W_lr
        center_pos = round((j-0.5)*sf);
        start = max(1, center_pos - kernel_radius);
        stop = min(W_original, center_pos + kernel_radius);
        
        kernel_start = max(1, kernel_radius + 2 - (center_pos - start));
        kernel_stop = min(length(v1), kernel_radius + 1 + (stop - center_pos));
        valid_kernel = v1(kernel_start:kernel_stop)';
        
        T(start:stop, j) = valid_kernel / sum(valid_kernel);
    end
end