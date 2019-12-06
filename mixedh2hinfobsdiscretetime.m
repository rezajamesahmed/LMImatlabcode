Ad = [1 2 3;
4 5 6;
7 8 9];
Bd1 = [0;0;1];
Bd2 = [0;1;1];
Cd2 = [1 0 0];
Cd1 = [1 0 1;
    1 1 1;
    0 2 3];
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
Mat1 = [P (P*Ad - Gd*Cd2) (P*Bd1 - Gd*Dd21);
    (P*Ad - Gd*Cd2)' P zeros(size((P*Bd1 - Gd*Dd21)));
    (P*Bd1 - Gd*Dd21)' zeros(size((P*Bd1 - Gd*Dd21)))' eye(1,1)]
Mat3 = [P (P*Ad - Gd*Cd2) (P*Bd1 - Gd*Dd21) zeros(3,1);
    (P*Ad - Gd*Cd2)' P zeros(3,1) Cd2';
    (P*Bd1 - Gd*Dd21)' zeros(1,3) gamma Dd21;
    zeros(1,3) Cd2 Dd21 gamma];
Mat2 = [Z P*Cd1;
    (P*Cd1)' P]

%define constraints
F = []
F = [F; Mat1>=smallnum*eye(size(Mat1))]
F = [F; Mat2>=smallnum*eye(size(Mat2))]
F = [F; Mat3>=smallnum*eye(size(Mat3))]
F = [F; trace(Z)<=v]
F = [F; P >=smallnum*eye(size(P))]
F = [F; Z >=smallnum*eye(size(Z))]

%optimize
solution = optimize(F,v+gamma)

%print optimal H2 observer
Ld = inv(value(P))*value(Gd)
H2norm = sqrt(value(v))
Hinfnorm = value(gamma)
