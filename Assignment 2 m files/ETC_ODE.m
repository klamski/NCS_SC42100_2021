function dx = ETC_ODE(t, x)
sigma = 0.5;
A = [5, 5.5; 0, -5];
B = [0; 1];
K = [10.1818 5];

xt = prevTrigState; %previous triggering state
e = xt - x;
dlmwrite('epsilon.txt', e, '-append');
if( norm(e, 1) > (1-sigma)*0.0131*norm(x,1) && t>0)
    setTrigState(x);
    event = 1;
    dlmwrite('ev_times.txt',t,'-append');
    dlmwrite('ev.txt',1,'-append')
else
    dlmwrite('ev.txt',0,'-append')
end

dx = (A-B*K)*x - B*K*e;