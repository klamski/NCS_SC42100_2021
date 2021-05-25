close all; clear; clc;
%% Ass2 Ex4  SubQ 1
A = [5, 5.5; 0, -5];
B = [0; 1];

S = [-11/20, 1; 1, 0];
J = [-5, 0; 0, 5];
S_inv = [0, 1; 1, 11/20];
Q_inv = S;
Q = S_inv;
% A = S*J*S_inv = Q_inv*J*Q;

syms h tau;
F0 = [exp(5*h) 11/20*(exp(5*h) - exp(-5*h)) 11/100*(exp(-5*h) + exp(5*h)); 0 exp(-5*h) -1/5*exp(-5*h); 0 0 0];
F1 = [0 0 -11/100; 0 0 1/5; 0 0 0];
F2 = [0 0 -11/100; 0 0 0; 0 0 0];

F = F0 + exp(-5*(h-tau))*F1 + exp(5*(h-tau))*F2;

G0 = [-22/100; 1/5; 1];
G1 = [11/100; -1/5; 0];
G2 = [11/100; 0; 0];

G = G0 + exp(-5*(h-tau))*G1 + exp(5*(h-tau))*G2;

%% Ex4 SubQ 2
K = [10.1818 5 0];
cnt = 1;
gamma = 0.01;
step = 0.01;
h_max = 1;
stable = zeros(1, h_max/step);
unstable = zeros(1, h_max/step);
tau_max = zeros(1, h_max/step);

for i = 0:step:h_max
h = i;
a1min = double(subs(exp(-5*h)));
a2min = double(subs(exp(5*h)));
a1max = 1;
a2max = 1;
%a1max = double(subs(exp(-5*0.9*h)));
%a2max = double(subs(exp(5*0.9*h)));
subsF0 = double(subs(F0));

% Stability analysis for different h values
P = sdpvar(3,3);
cons = [P >= 0, ... 
((subsF0 + a1min*F1 + a2min*F2)-(G0 + a1min*G1 + a2min*G2)*K)'*P*((subsF0 +a1min*F1 + a2min*F2)-(G0 + a1min*G1 + a2min*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1max*F1 + a2min*F2)-(G0 + a1max*G1 + a2min*G2)*K)'*P*((subsF0 + a1max*F1 + a2min*F2)-(G0 + a1max*G1 + a2min*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1min*F1 + a2max*F2)-(G0 + a1min*G1 + a2max*G2)*K)'*P*((subsF0 +a1min*F1 + a2max*F2)-(G0 + a1min*G1 + a2max*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1max*F1 + a2max*F2)-(G0 + a1max*G1 + a2max*G2)*K)'*P*((subsF0 + a1max*F1 + a2max*F2)-(G0 + a1max*G1 + a2max*G2)*K) - P + gamma*P <= 0];
obj = 0;
options = sdpsettings('verbose',0);
result = optimize(cons, obj, options);
disp(result.info)
Z = double(value(P));


try chol(Z)
    disp('Matrix is symmetric positive definite.')
    stable(cnt) = h;
    tau_max(cnt) = h;
catch ME
    disp('Matrix is not symmetric positive definite')
   % stable(cnt) = h;
   % tau_max(cnt) = find_max_tau(subsF0, F1, F2, G0, G1, G2, a1min, a2min, h, K, gamma);
end
    cnt = cnt + 1;
end

figure(1);
area(stable(1:end), tau_max(1:end))
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');

%% MORE LMIS
for i = 0:step:h_max
h = i;
a1min = double(subs(exp(-5*h)));
a2min = double(subs(exp(5*h)));
a1max = 1;
a2max = 1;
%a1max = double(subs(exp(-5*0.9*h)));
%a2max = double(subs(exp(5*0.9*h)));
subsF0 = double(subs(F0));

% Stability analysis for different h values
P = sdpvar(3,3);
cons = [P >= 0, ... 
((subsF0 + a1min*F1 + a2min*F2)-(G0 + a1min*G1 + a2min*G2)*K)'*P*((subsF0 +a1min*F1 + a2min*F2)-(G0 + a1min*G1 + a2min*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1max*F1 + a2min*F2)-(G0 + a1max*G1 + a2min*G2)*K)'*P*((subsF0 + a1max*F1 + a2min*F2)-(G0 + a1max*G1 + a2min*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1min*F1 + a2max*F2)-(G0 + a1min*G1 + a2max*G2)*K)'*P*((subsF0 +a1min*F1 + a2max*F2)-(G0 + a1min*G1 + a2max*G2)*K) - P + gamma*P <= 0, ...
((subsF0 + a1max*F1 + a2max*F2)-(G0 + a1max*G1 + a2max*G2)*K)'*P*((subsF0 + a1max*F1 + a2max*F2)-(G0 + a1max*G1 + a2max*G2)*K) - P + gamma*P <= 0];
obj = 0;
options = sdpsettings('verbose',0);
result = optimize(cons, obj, options);
disp(result.info)
Z = double(value(P));


try chol(Z)
    disp('Matrix is symmetric positive definite.')
    stable(cnt) = h;
    tau_max(cnt) = h;
catch ME
    disp('Matrix is not symmetric positive definite')
    stable(cnt) = h;
    tau_max(cnt) = find_max_tau(subsF0, F1, F2, G0, G1, G2, a1min, a2min, h, K, gamma);
end
    cnt = cnt + 1;
end
figure(2);
area(stable(1:end), tau_max(1:end))
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');