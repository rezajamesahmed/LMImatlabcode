%Discrete-Time H2-Optimal Observer Example, taken from (Caverly 5.1.2)
clear all
clc
smallnum = 0.00001;

%define matrices
Ad = [1 2 3;
    4 5 6;
    7 8 9]
Bd1 = [0 1;
    1 0;
    0 1]
Cd2 = [1 0 3;
    0 1 1]
Cd1 = [0 1 3;
    1 0 4;
    0 0 1]
Dd21 = [0 0;
    0 0]

%define decision variables
P = sdpvar(3,3)
Z = sdpvar(3,3)
Gd = sdpvar(3,2)
v = sdpvar(1)

%Build LMI
Mat1 = [P (P*Ad - Gd*Cd2) (P*Bd1 - Gd*Dd21);
    (P*Ad - Gd*Cd2)' P zeros(size((P*Bd1 - Gd*Dd21)));
    (P*Bd1 - Gd*Dd21)' zeros(size((P*Bd1 - Gd*Dd21)))' eye(2,2)]
Mat2 = [Z P*Cd1;
    (P*Cd1)' P]

%define constraints
F = []
F = [F; Mat1>=smallnum*eye(size(Mat1))]
F = [F; Mat2>=smallnum*eye(size(Mat2))]
F = [F; trace(Z)<=v]
F = [F; P >=smallnum*eye(size(P))]
F = [F; Z >=smallnum*eye(size(Z))]

%optimize
solution = optimize(F,v)

%print optimal H2 observer
Ld = inv(value(P))*value(Gd)
%print out H2 norm
mu = sqrt(value(v))

