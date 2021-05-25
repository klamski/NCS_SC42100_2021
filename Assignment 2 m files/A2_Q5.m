close all; clear; clc;
A = [5, 5.5; 0, -5];
B = [0; 1];
K = [10.1818 5];
P =  [ 356.9890  188.0931;  188.0931  124.4897];

state = [50; 100]; %initial state
setTrigState(state);
[t,xout] = ode45(@ETC_ODE, [0 5], state);


ev = dlmread('ev.txt');
times = dlmread('ev_times.txt');
%{
figure(1);
plot(t,xout)
xlabel("time");
ylabel("x state values");
%hold on

%figure(2);
%stem(times, ev)
%}
