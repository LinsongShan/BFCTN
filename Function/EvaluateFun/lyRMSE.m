function rmse = lyRMSE(Xfull, Xrecover)
% lyRMSE: Root mean Squared Error
%
% %[Syntax]%: 
%   rmse = lyRMSE(Xfull, Xrecover)
%  
% %[Inputs]%:
%   Xfull:       true data
%   Xrecover:    recover data
%
% %[Outputs]%:
%   rmse:        Root mean Squared Error
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  February 16, 2021
%

mse  = mean((Xfull(:) - Xrecover(:)).^2); % Mean Squared Error
rmse = sqrt(mse);                         % Root mean Squared Error