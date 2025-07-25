%% Function description
% The function performs backward substitution to solve Ux = y
% Inputs: U (an upper triangular matrix), y (a vector)
% Outputs: x (the solution vector)

%% Function code
function x = backwardSub(U,y)

% fetch the size of vector b
n = length(y);

% initialize the solution vector x
x = zeros(n,1);

% Solve for x using backward substitution
x(n) = y(n)/U(n,n);
for j=(n-1):-1:1
    x(j)=(y(j) -sum(U(j,j+1:end)'.*x(j+1:end)))/U(j,j);
end
