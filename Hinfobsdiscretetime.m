Ad = [1 2 3;
4 5 6;
7 8 9];
Bd1 = [0;0;1];
Cd1 = [1 0 0];
Dd21 = [3];
Dd22 = [0];
smallnum = 0.000001;

%define decision variables
P = sdpvar(3,3)
Z = sdpvar(3,3)
Gd = sdpvar(3,1)
v = sdpvar(1);
gamma = sdpvar(1);

%Build LMI
Mat3 = [P (P*Ad - Gd*Cd1) (P*Bd1 - Gd*Dd21) zeros(3,1);
    (P*Ad - Gd*Cd1)' P zeros(3,1) Cd1';
    (P*Bd1 - Gd*Dd21)' zeros(1,3) gamma Dd21;
    zeros(1,3) Cd1 Dd21 gamma];

%define constraints
F = []
F = [F; Mat3>=smallnum*eye(size(Mat3))]
F = [F; P >=smallnum*eye(size(P))]

%optimize
solution = optimize(F,gamma)

%print optimal H2 observer
Ld = inv(value(P))*value(Gd)
Hinfnorm = value(gamma)
