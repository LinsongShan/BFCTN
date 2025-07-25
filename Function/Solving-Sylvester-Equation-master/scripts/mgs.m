%% Function Description
% The function performs and returns the QR decomposition of a matrix
% It uses the Modified Gram-Schmidt procedure
% Input: Matrix A
% Output: Q and R (the QR decomposition of A

%% Function code
function [Q,R] = mgs(A)

% get the size of the matrix
n = size(A,1);

% initialize the Q and R matrices
Q = zeros(n,n);
R = zeros(n,n);

% store A in another matrix which would be changed
V=A;

% Apply Modified GS iteratively to compute Q and R
for j=1:n
    
    R(j,j) = norm(V(:,j));
    Q(:,j) = V(:,j)/R(j,j);
    
    for k=j+1:n
        R(j,k) = Q(:,j)'*V(:,k);
        V(:,k) = V(:,k)-R(j,k)*Q(:,j);
    end
    
end