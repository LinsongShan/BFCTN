function Y = Addnoise(X,SNR)
    DIM = size(X);                   
    sigma2 = var(X(:))*(1/(10^(SNR/10)));
    GN = sqrt(sigma2)*randn(DIM);
    Y = X + GN;
end