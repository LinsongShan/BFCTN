function HSI = generate_LRHSI_2(orgTensor,sf)
    s0 = sf/2;
    szX=size(orgTensor);
    % psf=fspecial('average',sf);
    psf        =    fspecial('gaussian',7,2);
    % psf = fspecial('motion', 4, 30);
    otf=psf2otf(psf,szX(1:2));
    H_otf=real(ifft2(fft2(orgTensor).*otf));
    HSI=H_otf(s0:sf:end,s0:sf:end,:);
end