function vpsnr = lyPSNR(Xfull, Xrecover)
% LYPSNR: Peak signal-to-noise ratio
%
% %[Syntax]%: 
%   vpsnr = lyPSNR(Xfull, Xrecover)
%  
% %[Inputs]%:
%   Xfull:       true data
%   Xrecover:    recover data
%
% %[Outputs]%:
%   vpsnr:       Peak signal-to-noise ratio
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  June 04, 2021
%

maxP = max(Xfull(:));

Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
[m,n,dim] = size(Xrecover);
MSE = norm(Xfull(:)-Xrecover(:))^2/(dim*m*n);
vpsnr = 10*log10(maxP^2/MSE);