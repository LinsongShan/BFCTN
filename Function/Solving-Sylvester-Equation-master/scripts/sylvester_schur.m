%% Program for solving Sylvester equation using Schur Factorization

% initialize the matrices for Sylvester equation
A =[ 9     2     3     0     1;
     2     7    -1     2    -1;
     3    -1     5    -1     3;
     0     2    -1     5     1;
     1    -1     3     1     9];
B =[ 13     1     5     4     3;
      1     5     0    -1     0;
      5     0     6     3     3;
      4    -1     3     5     4;
      3     0     3     4     9];
C =[21.6    1.3000   14.5000    6.1000    6.6000;
     1.4    6.3000   -2.2000    0.9000   -1.9000;
    15.6   -1.3000   17.7000    4.2000    7.2000;
     3.1    1.7000   -1.3000    5.6000    6.7000;
     8.8   -2.3000    8.3000    9.4000   24.8000];

% get the Schur factorization of A and B
[lambdaA,U] = QRapp1_modified(A, 30, 1);
[lambdaB,V] = QRapp1_modified(B, 30, 1);

% initialize Y
Y = zeros(5,5);

% set Z as U'CV
Z = U'*(C*V);

% solve for Y (U'XV) in lambdaA*Y + Y*lambdaB = Z
for i=1:5
    for j=1:5
        Y(i,j) = Z(i,j)/(lambdaA(i) + lambdaB(j));
    end
end

% solve for X in U'XV = Y, i.e. X = UYV'
X = U*Y*V';

% Solve the equation using in-built function
X_matlab = sylvester(A,B,C);

fprintf("\nThe error in solution using the developed method and inbuilt method is:\n\n");
disp(norm(X-X_matlab));
