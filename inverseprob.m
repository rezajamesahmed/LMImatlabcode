%define parameters
A=[1 2;
    0 2];
B = [0;
    1];
Q = [1 1;
    0 1];
R = [1 0; 
    0 1];

smallnum = 0.00001;

%Define variables
P = sdpvar(2,2);
P1 = sdpvar(2,2);
K = sdpvar(1,2);
gamma = sdpvar(1,1);

%define constraints
F=[];
F = [F; (A+B*K)'*P + P*(A+B*K) + K'*R*K + Q <= smallnum*eye(2,2)]
F = [F; B'*P+R*K <=smallnum*eye(2,2)]
F = [F; A'*P1+P1*A <= smallnum*eye(2,2)]

%optimize
solution = optimize(F,[])
