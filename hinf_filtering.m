%define parameters
A=[1 2 3;
    1 5 6;
    2 8 9];
B = [0; 1; 0];
C = [1 1 2];
D = [1];
L = [0 1 1];
smallnum = 0.00001;

%Define variables
R = sdpvar(3,3);
X = sdpvar(3,3);
M = sdpvar(3,3);
N = sdpvar(1,3);
Z = sdpvar(3,1);
Df = sdpvar(1,1);
gamma = sdpvar(1,1);

%define constraints
F=[];
F=[F;X>=smallnum*eye(size(X))];
F=[F;R-X>=smallnum*eye(size(X))];

Mat = [(R*A + A'*R + Z*C + C'*Z') (M'+Z*C+X*A)' (B'*R+D'*Z')' (L-Df*C)';
    (M'+Z*C+X*A) (M'+M) (B'*X+D'*Z')' (-N');
    (B'*R+D'*Z') (B'*X+D'*Z') (-gamma*eye(1,1)) (-Df*D)';
    (L-Df*C) -N (-Df*D) (-gamma*eye(1,1))];

F=[F;Mat<=-smallnum*eye(size(Mat))];

%optimize
solution = optimize(F,gamma)

%print Filter Coefficient Matrices
Af = inv(value(X))*value(M)
Bf = inv(value(X))*value(Z)
Cf = value(N)
Df = value(Df)

