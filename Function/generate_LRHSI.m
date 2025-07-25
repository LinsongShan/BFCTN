function HSI = generate_LRHSI(orgTensor, sf, kernel)
    % generate_LRHSI - Generate a low-resolution HSI with different blur kernels
    
    % kernel      - Degradation mode: 'gaussian', 'average', or 'motion'

    s0 = sf / 2;
    szX = size(orgTensor);

    % Select degradation kernel based on mode
    switch lower(kernel)
        case 'gaussian'
            psf = fspecial('gaussian', 4, 2);
        case 'average'
            psf = fspecial('average', sf);
        case 'motion'
            psf = fspecial('motion', 4, 30);
        otherwise
            error('Unknown degradation mode: %s. Choose ''gaussian'', ''average'', or ''motion''.', kernel);
    end

    otf = psf2otf(psf, szX(1:2));
    H_otf = real(ifft2(fft2(orgTensor) .* otf));
    HSI = H_otf(s0:sf:end, s0:sf:end, :);
end
