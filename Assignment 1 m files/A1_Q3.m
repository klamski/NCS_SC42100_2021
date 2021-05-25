close all; clear; clc;
A = [5, 5.5; 0, -5];
B = [0; 1];
C = [1 1];
D = 0;
if (rank(ctrb(A,B)) == rank(A))
    fprintf("The system is controllable \n");
end

%1) Pole placement
p = [-2 -3];
K = place(A, B, p)

syms h tau s;

Lambda = expm(A*h)
M0 = int(expm(A*s), s, 2*h - tau, h)*B
M1 = int(expm(A*s), s, 0, 2*h - tau)*B

F = [Lambda, M1, M0; [0 0], 0, 0; [0 0], 1, 0];
G = [0; 0; 1; 0];

Fcl = F-G*[K 0 0];

%% Stability Analysis for static controller 
step = 0.01;
cnt_c = 1;
cnt_r = 1;

for i=0:step:1
    for j=0:step:i
        h = i;
        tau = i+j;
        rho(cnt_r, cnt_c) = norm(max(abs(eig(double(subs(Fcl))))));
        cnt_c = cnt_c + 1; 
    end
    cnt_r = cnt_r + 1; 
    cnt_c = 1;
end

cnt_l = 1;
cnt_h = 1;
rhoT = rho';
first = 1;
cnt = 1;
for j = 1:size(rhoT,2)
    [row_l, col_l] = find((rhoT(:,j) > 0 & rhoT(:,j) <= 1),1,'first')
    [row_h, col_h] = find((rhoT(:,j) > 0 & rhoT(:,j) <= 1),1,'last')
    if(~isempty(row_l))
        hnew_l(cnt) = (j - 1)*step;
        taunew_l(cnt) = ((row_l - 1) + (j - 1))*step;
        hnew_h(cnt) = (j - 1)*step;
        taunew_h(cnt) = ((row_h - 1) + (j - 1))*step;
        cnt = cnt + 1;
    end
   
end

figure(2);
hold all
plot(hnew_l, taunew_l);
plot(hnew_h, taunew_h);
patch([hnew_l hnew_h], [taunew_l taunew_h], 'b')
hold off
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');

save('hnew_l.mat','hnew_l');
save('hnew_h.mat','hnew_h');
save('taunew_l.mat', 'taunew_l');
save('taunew_h.mat', 'taunew_h');
clear;

h_old_l = load('hnew_l.mat');
h_old_h = load('hnew_h.mat');
tau_old_l = load('taunew_l.mat');
tau_old_h = load('taunew_h.mat');
%{
h_old_l = hnew_l;
tau_old_l = taunew_l;
h_old_h = hnew_h;
tau_old_h = taunew_h;
%}
%% Stability analysis for dynamic controller 
% Exact Discrete Time System of Question 3 - tau in [h, 2h)
syms h;
F0big = [exp(5*h) 11/20*(exp(5*h) - exp(-5*h)) -22/100 -11/100*(exp(-5*h) - exp(5*h)); 0 exp(-5*h) 1/5 1/5*exp(-5*h); 0 0 0 0; 0 0 1 0];
F1big = [0 0 11/100*exp(-5*h) -11/100*exp(-5*h); 0 0 -1/5*exp(-5*h) 1/5*exp(-5*h); 0 0 0 0; 0 0 0 0];
F2big = [0 0 11/100*exp(5*h) -11/100*exp(5*h); 0 0 0 0; 0 0 0 0; 0 0 0 0];

Gbig = [0; 0; 1 ; 0];

step = 0.01;
h_max = 1;
cnt = 1;
gamma = 1;

for i=0:step:h_max
h = i;
sF0big = double(subs(F0big));
sF1big = double(subs(F1big));
sF2big = double(subs(F2big));

a13 = double(subs(exp(-5*0*h)));
a14 = double(subs(exp(5*0*h)));

a23 = double(subs(exp(5*h)));
a24 = double(subs(exp(-5*h)));

Y = sdpvar(4,4);    % Yj = Y for the case of  common quadratic lyapunov function
X = sdpvar(4,4);
Z = sdpvar(1,4);

cons = [Y>= 0, ... 
       %[X+X'-Y, X'*(sF0small + a11*F1small + a21*F2small)'-Z'*(G0small + a11*G1small + a21*G2small)'; (sF0small + a11*F1small + a21*F2small)*X-(G0small + a11*G1small + a21*G2small)*Z, gamma*Y]>= 0, ...
       %[X+X'-Y, X'*(sF0small + a12*F1small + a22*F2small)'-Z'*(G0small + a12*G1small + a22*G2small)'; (sF0small + a12*F1small + a22*F2small)*X-(G0small + a12*G1small + a22*G2small)*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*(sF0big + a13*sF1big + a23*sF2big)'-Z'*Gbig'; (sF0big + a13*sF1big + a23*sF2big)*X-Gbig*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*(sF0big + a13*sF1big + a24*sF2big)'-Z'*Gbig'; (sF0big + a13*sF1big + a24*sF2big)*X-Gbig*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*(sF0big + a14*sF1big + a23*sF2big)'-Z'*Gbig'; (sF0big + a14*sF1big + a23*sF2big)*X-Gbig*Z, gamma*Y]>= 0, ...
       [X+X'-Y, X'*(sF0big + a14*sF1big + a24*sF2big)'-Z'*Gbig'; (sF0big + a14*sF1big + a24*sF2big)*X-Gbig*Z, gamma*Y]>= 0 ];
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
    %eig1 = eig(-((sF0small + a11*F1small + a21*F2small)-(G0small + a11*G1small + a21*G2small)*Kc)'*P*((sF0small + a11*F1small + a21*F2small)-(G0small + a11*G1small + a21*G2small)*Kc) + P);
    %eig2 = eig(-((sF0small + a12*F1small + a22*F2small)-(G0small + a12*G1small + a22*G2small)*Kc)'*P*((sF0small + a12*F1small + a22*F2small)-(G0small + a12*G1small + a22*G2small)*Kc) + P);
    eig3 = eig(-((sF0big + a13*sF1big + a23*sF2big)-Gbig*Kc)'*P*((sF0big + a13*sF1big + a23*sF2big)-Gbig*Kc) + P);
    eig4 = eig(-((sF0big + a14*sF1big + a24*sF2big)-Gbig*Kc)'*P*((sF0big + a14*sF1big + a24*sF2big)-Gbig*Kc) + P);

    if(eig3 >= 0 & eig4 >= 0)%(eig1 >= 0 & eig2 >= 0 & eig3 >= 0 & eig4 >= 0)
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

%% Test new controller
syms h tau s;
A = [5, 5.5; 0, -5];
B = [0; 1];
Lambda = expm(A*h)
M0 = int(expm(A*s), s, 2*h - tau, h)*B
M1 = int(expm(A*s), s, 0, 2*h - tau)*B

F = [Lambda, M1, M0; [0 0], 0, 0; [0 0], 1, 0];
G = [0; 0; 1; 0];


Fcl = F-G*Kc;
step = 0.01;
cnt_c = 1;
cnt_r = 1;

for i=0:step:1
    for j=0:step:i
        h = i;
        tau = i+j;
        rho(cnt_r, cnt_c) = norm(max(abs(eig(double(subs(Fcl))))));
        cnt_c = cnt_c + 1; 
    end
    cnt_r = cnt_r + 1; 
    cnt_c = 1;
end

cnt_l = 1;
cnt_h = 1;
rhoT = rho';
first = 1;
cnt = 1;
for j = 1:size(rhoT,2)
    [row_l, col_l] = find((rhoT(:,j) > 0 & rhoT(:,j) <= 1),1,'first')
    [row_h, col_h] = find((rhoT(:,j) > 0 & rhoT(:,j) <= 1),1,'last')
    if(~isempty(row_l))
        hnew_l(cnt) = (j - 1)*step;
        taunew_l(cnt) = ((row_l - 1) + (j - 1))*step;
        hnew_h(cnt) = (j - 1)*step;
        taunew_h(cnt) = ((row_h - 1) + (j - 1))*step;
        cnt = cnt + 1;
    end
   
end

figure(3);
hold all
plot(hnew_l, taunew_l);
plot(hnew_h, taunew_h);
patch([hnew_l hnew_h], [taunew_l taunew_h], 'b')

plot(h_old_l.hnew_l, tau_old_l.taunew_l);
plot(h_old_h.hnew_h, tau_old_h.taunew_h);
patch([h_old_l.hnew_l h_old_h.hnew_h], [tau_old_l.taunew_l tau_old_h.taunew_h], 'g')

hold off
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');
legend('Controller 1', 'Controller 2');