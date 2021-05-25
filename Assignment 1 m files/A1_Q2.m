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

Fx = expm(A*h)
Fu = int(expm(A*s), s, h - tau, h)*B
G1 = int(expm(A*s), s, 0, h - tau)*B

F = [Fx Fu;0 0 0];
G = [G1; 1];

Fcl = F-G*[K 0]

rho = zeros(10, 10);
cnt_r = 1;
cnt_c = 1;
%{
for i=0:0.001:1
    for j = 0:0.001:i
        h = i;
        tau = j;
        rho(cnt_r, cnt_c) = norm(max(abs(eig(double(subs(Fcl))))));
        cnt_c = cnt_c + 1; 
    end
    cnt_r = cnt_r + 1;
    cnt_c = 1;
end

save('rho.mat', 'rho');

h_param(:) = 0:0.001:0.62;
tau_param(:) = 0:0.001:0.21;
ref = ones(63,22);
color = zeros(63, 22, 3);
color(:,:,1) = 255;
rr = rho(1:63, 1:22);
%{
figure(1);
surf(tau_param, h_param, rr); 
hold on;
surf(tau_param, h_param, ref, color);
ylabel('h');
xlabel('\tau');
zlabel('\rho(F_{cl}(h))');
%}
%[row, col] = find(0 < rho & rho <= 1);
%row = (row - 1)*0.01;
%col = (col - 1)*0.01;

rhoT = rho'
for j=1:size(rhoT,1)
for i = 1:j
if(rhoT(i,j) <= 1 && rhoT(i,j) ~= 0)
taunew(j) = (i - 1)*0.001;
hnew(j) = (j - 1)*0.001;
end
end
end

figure(2);
area(h_param(1,1:size(taunew,2)), taunew);
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');

h_no_delay = 0:0.001:((size(rho) - 1)*0.001);
ref_2d(1:size(h_no_delay,2)) = 1;
figure(3);
plot(h_no_delay, rho(:,1), 'LineWidth', 2);
hold on;
plot(h_no_delay, ref_2d,'r', 'Linewidth', 2);
xlabel('h');
ylabel('\rho(F_{cl}(h))');
title('Norm of the maximum eigenvalue for zero delays');
%}
%% New Controller
Fcl = F-G*[K 0.23]

rho_new_controller = zeros(10, 10);
cnt_r = 1;
cnt_c = 1;
step = 0.001;
h_param2(:) = 0:step:0.62;
tau_param2(:) = 0:step:0.21;
h_max = 0.45;

for i=0:step:h_max
    for j = 0:step:i
        h = i;
        tau = j;
        rho_new_controller(cnt_r, cnt_c) = norm(max(abs(eig(double(subs(Fcl))))));
        cnt_c = cnt_c + 1; 
    end
    cnt_r = cnt_r + 1;
    cnt_c = 1;
end

save('rho_new_controller.mat', 'rho_new_controller');

rhoT = rho_new_controller'
for j=1:size(rhoT,1)
for i = 1:j
if(rhoT(i,j) <= 1 && rhoT(i,j) ~= 0)
taunew2(j) = (i - 1)*step;
hnew2(j) = (j - 1)*step;
end
end
end

figure(4);
area(hnew2, taunew2);
xlabel('h');
ylabel('\tau');
title('Combinations of (h, \tau) retaining stability');
