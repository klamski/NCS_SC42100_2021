function tau = find_max_tau(subsF0, F1, F2, G0, G1, G2, a1min, a2min, h, K, gamma)

for tau = 0:0.001:h
a1max = double(subs(exp(-5*(h - tau))));
a2max = double(subs(exp(5*(h - tau))));

% Stability analysis for different tau in [0, h]
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
catch ME
    disp('Matrix is not symmetric positive definite')
    disp(tau);
    break;
end

end