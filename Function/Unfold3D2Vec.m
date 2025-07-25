function [Mat] = Unfold3D2Vec( X )

[I, J, K, L] = size(X);
Mat = zeros(L, I*J*K);

for l=1:L
    tmp = X(:,:,:,l);
    Mat(l,:) = tmp(:);
end