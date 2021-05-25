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
K = place(A, B, p);

syms h s;

F = expm(A*h);
G = int(expm(A*s), s, 0, h)*B;
Fcl = F-G*K
syms hl;
Fhl = expm(A*hl);
Ghl = expm(A*(hl - h))*int(expm(A*s), s, 0, h)*B;
Fclhl = Fhl - Ghl*K
epsilon = 0.001;
h = 0.1;
cnt = 0;

for hl = h:h:4*h    
    Dynamics = double(subs(Fclhl));
    P = sdpvar(2,2);
    Q = sdpvar(2,2);
    cons = [P>= epsilon, Q>= epsilon,...
            Dynamics'*P*Dynamics - P + Q <= 0];
    obj = 0;
    options = sdpsettings('verbose',0);
    result = optimize(cons, obj, options);
    disp(result.info);
    P = double(value(P));
    Q = double(value(Q));
    
    if(eig(P) >= 0 & eig(Q) >= 0)
        disp("P positive definite");
        cnt 
        cnt = cnt + 1;
    end    
end


disp("The maximum acceptable number of consecutive packet losses is:");
disp(cnt - 1);

%% Q3.2-3.3

step = 0.001;
h = 0.1;
A0 = double(subs(F-G*K));
A1 = double(subs(F));
cnt = 1;
epsilon = 0.00001;
for p = 0:step:1
    P = sdpvar(2,2);
    cons = [P>= epsilon, P-(1-p)*A0'*P*A0-p*A1'*P*A1 >= epsilon];
    obj = 0;
    options = sdpsettings('verbose',0);
    result = optimize(cons, obj, options);
    disp(result.info);
    P = double(value(P));
    if(eig(P) >= 0 & eig( P-(1-p)*A0'*P*A0-p*A1'*P*A1) > 0)
        disp("P positive definite");
        Pacc = P;
        p_acceptable(cnt) = p;
        Pbernoulli{cnt,1} = P;
        cnt = cnt + 1;        
    end 
end
index = find(p_acceptable == max(p_acceptable(:)));
disp("the maximimum probability for Bernoulli MSS");
p = p_acceptable(index)
disp("for the matrix P ");
P = Pbernoulli{index,1}

%% Q 3.4

step = 0.001;
h = 0.1;
F0cl = double(subs(F-G*K));
F1cl = double(subs(F));
epsilon = 0.001;
cnt = 1;
for lamda = 0:step:1
    P = sdpvar(2,2);
    cons = [P >= epsilon, -(F0cl'*P*F0cl - lamda*P)>= 0];
    obj = 0;
    options = sdpsettings('verbose',0);
    result = optimize(cons, obj, options);
    disp(result.info);
    P = double(value(P));
    
    if(eig(P)>0 & eig(-(F0cl'*P*F0cl - lamda*P)) >= 0)
        L = sdpvar(1,1);
        cons = [L >= epsilon,  -(F1cl'*P*F1cl - L*P) >= 0];
        obj = 0;
        options = sdpsettings('verbose',0);
        result = optimize(cons, obj, options);
        disp(result.info);
        L = double(value(L));
        if(eig(-(F1cl'*P*F1cl - L*P))>=0)
            for p = 0:step:1
                if((lamda^(1-p))*(L^p)<=1)                    
                    ASS_star(cnt,1) = p;
                    ASS_star(cnt,2) = lamda;
                    ASS_star(cnt,3) = L;
                    Pmatrices{cnt,1} = P;
                    cnt = cnt + 1;
                    ASS_star
                end
            end
        end
    end
end

index = find(ASS_star(:,1) == max(ASS_star(:,1)));
disp("upper bound p star ASS equals");
p = ASS_star(index, 1)
disp("for lambda");
lamda = ASS_star(index, 2)
disp("for L");
L = ASS_star(index,3)
disp("and P matrix for CQLF");
P = Pmatrices{index, 1}

%% Q 3.5
step = 0.01;
h = 0.1;
A0 = double(subs(F-G*K));
A1 = double(subs(F));
cnt = 1;
cnt_c = 2;
epsilon = 0.1;


for p00 = 0.7:step:1    
    %p00 = sdpvar(1,1);
    %p11 = sdpvar(1,1);
    for p11 = 0:step:1
        P0 = sdpvar(2,2);
        P1 = sdpvar(2,2);
        cons = [P0>= epsilon, P1>= epsilon, ... %p11 >= 0,  p11 <= 1, ... , p00 >= epsilon,p00 <= 1 ...
                P0-p00*A0'*P0*A0-(1-p00)*A1'*P1*A1 >= epsilon, ...
                P1-(1-p11)*A0'*P0*A0-p11*A1'*P1*A1 >= epsilon];
        obj = 0;
        options = sdpsettings('verbose',0);
        result = optimize(cons, obj, options);
        disp(result.info);
        %p11 = value(p11)
        P0 = value(P0)
        P1 = value(P1)
        if(eig(P0) > 0 & eig(P1) > 0 & eig( P0-p00*A0'*P0*A0-(1-p00)*A1'*P1*A1 )>0 & eig(P1-(1-p11)*A0'*P0*A0-p11*A1'*P1*A1)>0 )            
                pairs(cnt, 1) = p00;
                pairs(cnt, cnt_c) = p11;
                cnt_c = cnt_c+ 1;           
        end
    end
    cnt = cnt + 1;
    cnt_c = 2;
end

figure(1);
scatter(pairs(:,1), pairs(:,2:end), 'k');
xlabel('p_{00}');
ylabel('p_{11}');
title('Combinations of p00 and p11 retaining MSS');
grid;

%{
for i = 1:size(pairs,1)
    p11max(i,1) = max(pairs(i, 2:end));
end
pairs(1:10,1) = 0.8;
figure(2);
area(pairs(:,1), p11max(:,1));
xlabel('p_{00}');
ylabel('p_{11}');
title('Combinations of p00 and p11 retaining MSS');
grid;
%}


