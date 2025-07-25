%% Function description
% The function reduces a matrix to upper Hessenburg form
% using Householder reflections
%
% Inputs: A (matrix)
%
% Outputs: H (the upper Hessenberg form of A)
%          T (the transformation matrix)

%% Function code
function [H,T] = hessenberg_householder(A)

% get the dimensions of A
m = length(A);

% initialize the transformation matrix
T = eye(m);

H=A;
 
for k = 1:m-2
    
    vk = H(k+1:m,k);
    alpha = sign(vk(1))*norm(vk);
    vk(1) = vk(1) - alpha;
    vk = vk/norm(vk);
    H(k+1:m,k:m) = H(k+1:m,k:m) - 2*vk*(vk.'*H(k+1:m,k:m));
    H(1:m,k+1:m) = H(1:m,k+1:m) - 2*(H(1:m,k+1:m)*vk)*vk.';
    Q = eye(m);
    Q(k+1:m,k+1:m) = Q(k+1:m,k+1:m) - 2*(vk*vk');
    T = T*Q;

end