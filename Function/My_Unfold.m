function [Mat] = My_Unfold( X )

[I, J, K, L] = size(X);
Mat = zeros(L, I*J*K);

for l=1:L
    tmp = X(:,:,:,l);
    tmp1 = permute(tmp,[3,2,1]);
    Mat(l,:) = tmp1(:);
end