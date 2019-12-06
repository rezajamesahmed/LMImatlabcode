%define params
A0=[1 2 3;
    4 5 6;
    7 8 9];
A1=[3 8 4;
    5 6 9;
    1 3 7];
A2=[0 1 3;
    1 1 1;
    4 7 9];

delA0=[.01 0 0;
    .02 0 .03;
    0 0 .009];
delA1=[.02 .02 .01;
    .02 0 .01;
    0 .02 .01];
delA2=[0 .01 .03;
    .0013 0 .01;
    0.3 .07 .09];

B =[0 1; 0 10; 3 0]
Cd= [1 2 5;
    0 0 1];
Cp=[50 1 6;
    0 2 0];
smallnum = 0.000001;

%define variables
Kp = sdpvar(2,2);
Kd = sdpvar(2,2);
eps0 = det(A0);
eps1 = det(A1);

%define constraints
F = [];

F = [F; A0 - B*Kd*Cd >= eps0*eye(size(A0))];
F = [F; (A1 - B*Kp*Cp)+(A1-B*Kp*Cp)' >=eps1*eye(size(A0))];

%optimize
solution = optimize(F,[])

%print Kp and Kd
Kpmat = value(Kp)
Kdmat = value(Kd)