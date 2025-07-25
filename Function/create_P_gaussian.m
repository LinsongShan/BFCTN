function S  = create_P_gaussian(patchsize, space_ratio, psf, psf_size)
        [U, S, V] = svd(psf);
        u1 = U(:, 1); % 高度方向滤波向量
        v1 = V(:, 1)'; % 宽度方向滤波向量
        
        % 归一化（保证能量守恒）
        u1 = u1 / sum(u1);
        v1 = v1 / sum(v1);
        
        % 构造降采样矩阵 S (高度方向)
        H_original = patchsize; % 原始图像高度
        H_lr = H_original / space_ratio;
        S = zeros(H_lr, H_original);
        for i = 1:H_lr
            center = (i-1)*space_ratio + (space_ratio/2); % 下采样中心位置
            start = max(1, round(center - (psf_size-1)/2));
            stop = min(H_original, round(center + (psf_size-1)/2));
            S(i, start:stop) = u1(1:(stop-start+1));
        end
        
        % % 构造降采样矩阵 T (宽度方向)
        % W_original = patchsize; % 原始图像宽度
        % W_lr = W_original / space_ratio;
        % T = zeros(W_original, W_lr);
        % for j = 1:W_lr
        %     center = (j-1)*space_ratio + (space_ratio/2); % 下采样中心位置
        %     start = max(1, round(center - (psf_size-1)/2));
        %     stop = min(W_original, round(center + (psf_size-1)/2));
        %     T(start:stop, j) = v1(1:(stop-start+1))';
        % end
end