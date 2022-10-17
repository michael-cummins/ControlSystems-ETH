close all
clear all
clc

% define constants
z_d = 2;
psi_d = pi/4;
m = 0.03;
g = 9.81;
Ix = 1.5e-5;
Iy = Ix;
Iz = 3e-5;
kx = 4.5e-3;
ky = 4.5e-3;
kz = 4.5e-3;
kp = 4.5e-4;
kq = 4.5e-4;
kr = 4.5e-4;

% Define state space model
A = [0 1; 0 -15];
B = [0; 3.3333e+04];
C = [1 0];
D = 0;
T = 0.01;

% Discretize the system
[Ad, Bd, Cd, Dd] = c2dm(A,B,C,D,T,'zoh');
sys = ss(A, B, C, D);

% Initial Conditions
x0 = [0;0];
t = (0:0.01:1);

% impulse response evolution
u4 = (t==0)*3e-05;
y1 = lsim(sys, u4, t, x0);

% step response evolution
u4 = (t>=0)*3e-05;
y2 = lsim(sys, u4, t, x0);

% sine wave response evolution
u4 = sin(t)*3e-05;
y3 = lsim(sys, u4, t, x0);

% Plotting results















