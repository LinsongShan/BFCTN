function srr = lySuccRecR(Xfull, Xrecover, tau)
% lySuccRecoR: Successful Recovery Ratio
%
% %[Syntax]%: 
%   srr = lySuccRecR(Xfull, Xrecover)
%  
% %[Inputs]%:
%   Xfull:       true data
%   Xrecover:    recover data
%   tau:         threshold value
%
% %[Outputs]%:
%   srr:        Successful Recovery Ratio
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  February 16, 2021
%

% ref: Xie K, Wang L, Wang X, et al. Accurate recovery of internet traffic data: A sequential tensor completion approach. 
%      IEEE/ACM Transactions on Networking, 2018, 26(2): 793-806.


absDif = abs(Xfull-Xrecover);
% absDif = abs(Xfull-Xrecover)./Xrecover;
O = absDif < tau;
srr = sum(O(:))/numel(O);

