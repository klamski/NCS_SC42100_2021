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

%2)exact discrete time model
syms h s;
Atil = [A B; 0 0 0]
a = expm(Atil*h)
F = [a(1,1) a(1,2); a(2,1) a(2,2)];
G = [a(1,3); a(2,3)];
Fcl = F-G*K

Fcln = expm(A*h) - int(expm(A*s),s, 0, h)*B*K
%3) Stability analysis as a function of the sampling interval h

j = 1;
for i = 0:0.01:10
    h = i;
    x_ax(j) = h;
    Fcl2 = Fcl;
    rho(j) = norm(max(abs(eig(double(subs(Fcl2))))));    
    j = j + 1;
end

ref(1:size(rho,2)) = 1;
stop = 70;%size(rho,2);
%plot
figure(1);
plot(x_ax(1,1:stop), rho(1,1:stop), 'b', 'LineWidth', 2);
hold on;
plot(x_ax(1,1:stop), ref(1,1:stop), 'r', 'LineWidth', 2);
xlabel('h');
ylabel('\rho(F_{cl}(h))');
title('Norm of the maximum eigenvalue');

%[row, col] = find(rho == 1)
