close all; clear; clc;
A = [5, 5.5; 0, -5];
B = [0; 1];
C = [1 1];
D = 0;

%% Exact Discrete Time System of Question 2 - tau in [0, h)
syms h tau s;

Fx = expm(A*h);
Fu = int(expm(A*s), s, h - tau, h)*B;
G1 = int(expm(A*s), s, 0, h - tau)*B;

Fsmall = [Fx Fu;0 0 0];
Gsmall = [G1; 1];

%% Exact Discrete Time System of Question 3 - tau in [h, 2h)

Lambda = expm(A*h);
M0 = int(expm(A*s), s, 2*h - tau, h)*B;
M1 = int(expm(A*s), s, 0, 2*h - tau)*B;

Fbig = [Lambda, M1, M0; [0 0], 0, 0; [0 0], 1, 0];
Gbig = [0; 0; 1; 0];

%% extended state vector xe(k) = [x1(k) x2(k) u(k-1) u(k-2)]
Fsmall = [Fx Fu [0;0];0 0 0 0; 0 0 1 0];
Gsmall = [G1; 1; 0];

%% Controller Synthesis of the Switched Model - tau in [0, 2h)
tau = 0.2*h;
F1small = subs(Fsmall);
G1small = subs(Gsmall);

tau = 0.5*h;
F2small = subs(Fsmall);
G2small = subs(Gsmall);

tau = h;
F1big = subs(Fbig);
G1big = subs(Gbig);

tau = 1.5*h;
F2big = subs(Fbig);
G2big = subs(Gbig);


% LMIS Theorem 7.2 from Book
step = 0.01;
h_max = 1;
cnt = 1;
gamma = 1;
for i=0:step:h_max
h = i;
sF1small = double(subs(F1small)); sG1small = double(subs(G1small));
sF2small = double(subs(F2small)); sG2small = double(subs(G2small));
sF1big = double(subs(F1big)); sG1big = double(subs(G1big));
sF2big = double(subs(F2big)); sG2big = double(subs(G2big));

Y = sdpvar(4,4);    % Yj = Y for the case of  common quadratic lyapunov function
X = sdpvar(4,4);
Z = sdpvar(1,4);

cons = [Y>= 0, ... 
       [X+X'-Y, X'*sF1small'-Z'*sG1small'; sF1small*X-sG1small*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*sF2small'-Z'*sG2small'; sF2small*X-sG2small*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*sF1big'-Z'*sG1big'; sF1big*X-sG1big*Z, gamma*Y]>= 0 ...%, ...
       %[X+X'-Y, X'*sF2big'-Z'*sG2big'; sF2big*X-sG2big*Z, gamma*Y]>= 0 ...
       ];
obj = 0;
options = sdpsettings('verbose',0);
result = optimize(cons, obj, options);
disp(result.info)
Y = (value(Y));
if(eig(Y) >= 0)
    stable(cnt) = h;
    Ksynth(cnt,1:4) = double(Z)*inv(double(X));
    Kc = Ksynth(cnt,1:4);
    disp('Matrix is symmetric positive definite.')
    
    P = inv(double(Y));
    eig1 = eig(-(sF1small - sG1small*Kc)'*P*(sF1small - sG1small*Kc) + P);
    eig2 = eig(-(sF2small - sG2small*Kc)'*P*(sF2small - sG2small*Kc) + P);
    eig3 = eig(-(sF1big - sG1big*Kc)'*P*(sF1big - sG1big*Kc) + P);
    %eig4 = eig(-(sF2big - sG2big*Kc)'*P*(sF2big - sG2big*Kc) + P);

    if(eig1 >= 0 & eig2 >= 0 & eig3 >= 0 );%& eig4 >= 0)
        disp("OK");
        stability(cnt) = h;
    end
end

cnt = cnt + 1;   
end

fprintf("\n The controller which maximizes the sampling interval h is: \n")
Kc = Ksynth(find(stability == max(stability)), : )
fprintf("\n For the following values of sampling interval h: \n");
stability


%% Verification of Results
%{
Kc = Ksynth(14,:);
step = 0.01;
h_max = 1;
cnt = 1;
for i = 0:step:h_max
    h = i;
    sF1small = double(subs(F1small)); sG1small = double(subs(G1small));
    sF2small = double(subs(F2small)); sG2small = double(subs(G2small));
    sF1big = double(subs(F1big)); sG1big = double(subs(G1big));
    sF2big = double(subs(F2big)); sG2big = double(subs(G2big));
    
    P = sdpvar(4,4);
    Q = sdpvar(4,4);
    cons = [ P >= 0, Q >= 0, [(sF1small - sG1small*Kc)'*P*(sF1small - sG1small*Kc) - P + Q ]<=0, [(sF2small - sG2small*Kc)'*P*(sF2small - sG2small*Kc) - P + Q] <= 0, [(sF1big - sG1big*Kc)'*P*(sF1big - sG1big*Kc) - P + Q] <= 0];%, [(sF2big - sG2big*Kc)'*P*(sF2big - sG2big*Kc) - P + Q] <= 0 ];
    obj = 0;
    options = sdpsettings('verbose',0);
    result = optimize(cons);%, obj, options);
    disp(result.info)
    P = (value(P))  
    
    if(eig(P) >= 0 )
        disp("OK");
        stability2(cnt) = h;
    end
    cnt = cnt + 1;
end
%}