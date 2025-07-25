function [S, T] = create_P_motion(patchsize, sf, psf)
    % create_P_motion - Approximate motion blur + downsampling operator
    % as separated height and width matrices S and T.
    %
    % Inputs:
    %   patchsize : size of original patch (assumed square)
    %   sf        : downsampling factor
    %   psf       : motion blur kernel (e.g. fspecial('motion', len, theta))
    %
    % Outputs:
    %   S : downsampling matrix for height dimension (size floor(patchsize/sf) × patchsize)
    %   T : downsampling matrix for width dimension (size patchsize × floor(patchsize/sf))
    
    % --- Step 1: Symmetrize psf to square matrix ---
    [h, w] = size(psf);
    if h == 1 || w == 1
        % vector psf, embed into square matrix preserving orientation
        kernel_length = max(h, w);
        symmetric_psf = zeros(kernel_length, kernel_length);
        center = ceil(kernel_length / 2);
        if h > w
            symmetric_psf(:, center) = psf(:);
        else
            symmetric_psf(center, :) = psf(:)';
        end
    else
        % non-square psf: resize to square via bilinear interp
        [X, Y] = meshgrid(1:w, 1:h);
        new_size = max(h, w);
        [XI, YI] = meshgrid(linspace(1, w, new_size), linspace(1, h, new_size));
        symmetric_psf = interp2(X, Y, psf, XI, YI, 'linear', 0);
    end
    
    % --- Step 2: SVD and retain main components ---
    [U, S_mat, V] = svd(symmetric_psf);
    singular_vals = diag(S_mat);
    energy_ratio = singular_vals(1) / sum(singular_vals);
    
    if energy_ratio < 0.85 && length(singular_vals) > 1
        % Use top two components for better approx
        u1 = U(:,1) * sqrt(singular_vals(1)) + U(:,2) * sqrt(singular_vals(2));
        v1 = V(:,1)' * sqrt(singular_vals(1)) + V(:,2)' * sqrt(singular_vals(2));
    else
        % Use rank-1 approx
        u1 = U(:,1) * sqrt(singular_vals(1));
        v1 = V(:,1)' * sqrt(singular_vals(1));
    end
    
    % Normalize kernels
    u1 = u1 / sum(u1);
    v1 = v1 / sum(v1);
    
    % --- Step 3: Construct downsampling matrices with corrected indexing ---
    H_orig = patchsize;
    H_lr = floor(H_orig / sf);
    S = zeros(H_lr, H_orig);  % height downsampling
    
    % For u1 (height kernel)
    center_kernel_u = floor(length(u1)/2) + 1;   % Center index of kernel
    left_radius_u = center_kernel_u - 1;          % Left radius (number of elements left of center)
    right_radius_u = length(u1) - center_kernel_u; % Right radius (number of elements right of center)
    
    for i = 1:H_lr
        center_pos = round((i-0.5) * sf);
        % Calculate spatial bounds using kernel's actual radii
        start_pos = max(1, center_pos - left_radius_u);
        stop_pos = min(H_orig, center_pos + right_radius_u);
        
        left_extent = center_pos - start_pos;   % Actual left extent in spatial domain
        right_extent = stop_pos - center_pos;   % Actual right extent in spatial domain
        
        % Calculate kernel segment indices
        kernel_start = center_kernel_u - left_extent;
        kernel_stop = center_kernel_u + right_extent;
        
        valid_kernel = u1(kernel_start:kernel_stop);
        % Assign normalized kernel (ensure sizes match)
        S(i, start_pos:stop_pos) = valid_kernel' / sum(valid_kernel);
    end
    
    W_orig = patchsize;
    W_lr = floor(W_orig / sf);
    T = zeros(W_orig, W_lr);  % width downsampling
    
    % For v1 (width kernel)
    center_kernel_v = floor(length(v1)/2) + 1;   % Center index of kernel
    left_radius_v = center_kernel_v - 1;          % Left radius
    right_radius_v = length(v1) - center_kernel_v; % Right radius
    
    for j = 1:W_lr
        center_pos = round((j-0.5) * sf);
        % Calculate spatial bounds using kernel's actual radii
        start_pos = max(1, center_pos - left_radius_v);
        stop_pos = min(W_orig, center_pos + right_radius_v);
        
        left_extent = center_pos - start_pos;
        right_extent = stop_pos - center_pos;
        
        % Calculate kernel segment indices
        kernel_start = center_kernel_v - left_extent;
        kernel_stop = center_kernel_v + right_extent;
        
        valid_kernel = v1(kernel_start:kernel_stop)';
        % Assign normalized kernel (ensure sizes match)
        T(start_pos:stop_pos, j) = valid_kernel / sum(valid_kernel);
    end
end