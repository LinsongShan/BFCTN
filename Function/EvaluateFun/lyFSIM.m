function FSIM = lyFSIM(imagery1, imagery2)
%==========================================================================
%
% Input:
%   imagery1 - the reference tensor
%   imagery2 - the target tensor

%
% Output:
%   fsim - Feature SIMilarity
%==========================================================================


[~,~,l] = size(imagery1);
fsim = zeros(1, l);

for i = 1:l 
    fsim(i) = FeatureSIM(imagery1(:, :, i), imagery2(:, :, i));
end

FSIM = mean(fsim);

end