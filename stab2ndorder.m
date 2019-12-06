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
B =[0; 10; 0]
Cd= [1 2 5];
Cp=[50 1 6];
smallnum = 0.000001;

%define variables
Kp = sdpvar(3,3);
Kd = sdpvar(3,3);

%define constraints
F = [];
F = [F; A2>0;A1+A1'>0;A0>0];
F = [F; A0 - B*Cd*Kd >= smallnum*eye(size(A0))];
F = [F; (A1 - B*Cp*Kp)+(A1-B*Cp*Kp)' >=smallnum*eye(size(A0))];

%optimize
solution = optimize(F,[])

%print Kp and Kd
Kpmat = value(Kp)
Kdmat = value(Kd)