function rrse = lyRRSE(Xfull, Xrecover)
% lyRRSE: Root Relative Square Error
%
% %[Syntax]%: 
%   rrse = lyRRSE(Xfull, Xrecover)
%  
% %[Inputs]%:
%   Xfull:       true data
%   Xrecover:    recover data
%
% %[Outputs]%:
%   rrse:        Root Relative Square Error
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  February 16, 2021
%


rrse = sqrt(sum((Xrecover(:) - Xfull(:)).^2)/sum(Xfull(:).^2)); % Root Relative Square Error