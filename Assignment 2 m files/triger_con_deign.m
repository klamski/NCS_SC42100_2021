close all; clear; clc;
A = [5, 5.5; 0, -5];
B = [0; 1];
K = [10.1818 5];

P = sdpvar(2,2);
Q = sdpvar(2,2);
epsilon = 0.001;
cons = [P >= epsilon*eye(2), Q >= epsilon*eye(2,2), ...
        ((A-B*K)'*P+P*(A-B*K)) + Q == 0];
obj = 0;
result = optimize(cons, obj);
disp(result.info);
P = value(P)
Q = value(Q)

min(eig(Q))/(2*norm(P*B*K,1))