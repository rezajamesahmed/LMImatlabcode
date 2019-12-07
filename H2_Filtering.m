%define parameters
A=[1 2 3;
    0 2 6;
    2 7 3];
B = [0 1; 1 2; 0 2];
C = [1 1 2;
    0 0 1];
D = [1 0; 2 1];
L = [0 1 1;
    1 3 2];
smallnum = 0.00001;

%Define variables
R = sdpvar(3,3);
X = sdpvar(3,3);
M = sdpvar(3,3);
N = sdpvar(2,3);
Z = sdpvar(3,2);
Q = sdpvar(2,2);
gamma = sdpvar(1,1);

%define constraints
F=[];
F=[F;R-X>=smallnum*eye(size(X))];
F=[F;Q<=gamma];
Mat = [(R*A + A'*R + Z*C + C'*Z') (M'+Z*C+X*A)' (B'*R+D'*Z')';
    (M'+Z*C+X*A) (M'+M) (B'*X+D'*Z')' ;
    (B'*R+D'*Z') (B'*X+D'*Z') (-eye(2,2))];

Mat2 = [-Q L -N;
    L' -R -X';
    -N' -X -X]

F=[F;Mat<=-smallnum*eye(size(Mat))];
F=[F;Mat2<=-smallnum*eye(size(Mat2))];

%optimize
solution = optimize(F,gamma)

%print preliminary Filter Coefficient Matrices
Af = inv(value(X))*value(M)
Bf = inv(value(X))*value(Z)
Cf = value(N)

%print completed Filter Coefficient Matrices
Atilde = [A zeros(size(A));
    Bf*C Af]
Btilde = [B;Bf*D]
Ctilde = [L -Cf]