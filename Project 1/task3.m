close all
clear all
clc

% define constants
z_d = 2;
psi_d = pi/4;
m = 0.03;
g = 9.81;
Ix = 1.5e-5;
Iy = 1.5e-5;
Iz = 3e-05;
kx = 4.5e-03;
ky = 4.5e-03;
kz = 4.5e-03;
kp = 4.5e-04;
kq = 4.5e-04;
kr = 4.5e-04;

% Define state space model
A = [0 1; 0 -15];
B = [0; 3.3333e+04];
C = [1 0];
D = 0;
T = 0.01;

% Discretize the system
% [Ad, Bd, Cd, Dd] = c2dm(A,B,C,D,T,'zoh');
sys = ss(A, B, C, D);

% Initial Conditions
x0 = [0;0];
t = (0:T:1);

% impulse response evolution
y1 = impulse(sys, t)*3e-05;


% step response evolution
y2 = step(sys, t)*3e-05;

% sine wave response evolution
u4 = sin(t)*3e-05;
y3 = lsim(sys, u4, t, x0);

% Analysing transfer function
Hs = tf(sys);


% Plotting results
subplot(3,1,1);
plot(t, y1)
title('Impulse Response')
xlabel('Time')
ylabel('Output')
axis([0 1 0 0.08])
subplot(3,1,2)
plot(t, y2)
title('Step Response')
xlabel('Time')
ylabel('Output')
axis([0 1 0 0.08])
subplot(3,1,3)
plot(t, y3)
title('Sine wave repsonse')
axis([0 1 0 0.03])
xlabel('Time')
ylabel('Output')













