function mae = lyMAE(Xfull, Xrecover)
% lyMAE: Average of the Absolute Error
%
% %[Syntax]%: 
%   mae = lyMAE(Xfull, Xrecover)
%  
% %[Inputs]%:
%   Xfull:       true data
%   Xrecover:    recover data
%
% %[Outputs]%:
%   mae:        Average of the Absolute Error
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  February 16, 2021
%

% ref: Xie K, Wang L, Wang X, et al. Accurate recovery of internet traffic data: A sequential tensor completion approach. 
%      IEEE/ACM Transactions on Networking, 2018, 26(2): 793-806.

absDif = abs(Xfull-Xrecover);
mae = sum(absDif(:))/numel(Xrecover);