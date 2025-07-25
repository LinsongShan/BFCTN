function S  = create_P_gaussian(patchsize, space_ratio, psf, psf_size)
        [U, S, V] = svd(psf);
        u1 = U(:, 1); 
        v1 = V(:, 1)'; 

        u1 = u1 / sum(u1);
        v1 = v1 / sum(v1);
        
        H_original = patchsize; 
        H_lr = H_original / space_ratio;
        S = zeros(H_lr, H_original);
        for i = 1:H_lr
            center = (i-1)*space_ratio + (space_ratio/2); % 下采样中心位置
            start = max(1, round(center - (psf_size-1)/2));
            stop = min(H_original, round(center + (psf_size-1)/2));
            S(i, start:stop) = u1(1:(stop-start+1));
        end

end