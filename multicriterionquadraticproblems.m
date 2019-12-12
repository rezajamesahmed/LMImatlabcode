%Multicreiterion LQG

A = [1 2 3;
    4 5 7;
    7 8 9];
B = [1 2 3;
    3 4 3;
    5 6 3]
C = [1 0 1;
    0 1 1;
    3 3 3];
D = [0 0 0 ;
    0 0 0;
    0 0 0];

Q = sdpvar(3,3);
R = sdpvar(3,3);
C = sdpvar(3,3);
S = sdpvar(3,3);
P = sdpvar(3,3);
gamma = sdpvar(1,1);
smallnum = 0.000001;



for i= 1:1:3
Q = Q + sdpvar(3,3)*i;
R = R + sdpvar(3,3)*i;
C = C + sdpvar(3,3)*i;
S = S + sdpvar(3,3)*i;
gamma = sdpvar(1,1);
  
end



F = [];
F=[F;[(A'*P + P*A + Q) (P*Q+C');(B'*P+C) (R)] >=smallnum];


%optimize
solution = optimize(F,gamma)



