%% Function description
% The function solves the Sylvester equation
% using Kronecker product and QR algorithm
%
% Inputs: A, B, C (matrices of the Sylvester equation)
%
% Outputs: X (the solution to the equation)

%% Function code
function X = sylvester_kronecker(A, B, C)

% get the dimensions of A and B
dimA = length(A);
dimB = length(B);

% create a linear system using Kronecker product and vectorization
K = kron(eye(dimA),A) + kron(B',eye(dimB));
vecC = C(:);

% get the QR decomposition of K using mgs
[Q,R] = mgs(K);

% get the system in form of Rx = b, where b = Q'*vecC
b = Q'*vecC;

%solve Rx = b using backward substitution
vecX = backwardSub(R,b);

X = reshape(vecX, dimA, dimB);



 