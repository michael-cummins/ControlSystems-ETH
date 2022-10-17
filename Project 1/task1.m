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

% equilibrium points
sys1_eq = [0 0];
sys2_eq = [0 0 0 0 0 psi_d m*g];
sys3_eq = [0 0 0 0 0 0];

% system 1 linearization
syms phi theta 
sys1_states = [phi theta];
sys1_A = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
          0 cos(phi) -sin(phi);
          0 sin(phi)/cos(theta) cos(phi)*cos(theta)];
sys1_A_lin = subs(sys1_A, sys1_states, sys1_eq);

% system 2 linearization
syms xdot ydot zdot phi theta psi u1

states = [xdot ydot zdot];
euler_angles = [phi theta psi]; 
input = u1;

xddot = (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*u1 - kx*xdot;
yddot = (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*u1 - ky*ydot;
zddot = (cos(phi)*cos(theta)*u1 - m*g - kz*zdot);
sys2 = [xddot; yddot; zddot]/m;

sys2_Alin = subs(jacobian(sys2, states), [states euler_angles input], sys2_eq);
sys2_Blin = subs(jacobian(sys2, input), [states euler_angles input], sys2_eq);

A2 = double(vpa(sys2_Alin));
B2 = double(vpa(sys2_Blin));

% system 3 linearization
syms p q r u2 u3 u4
states = [p q r];
inputs = [u2 u3 u4];

pdot = ((Ix - Iz)*q*r + u2 -kp*p)/Ix;
qdot = ((Iz - Ix)*p*r + u3 - kq*q )/Iy;
rdot = ((Ix - Iy)*p*q + u4 - kr*r )/Iz;
sys3 = [pdot; qdot; rdot];

sys3_Alin = subs(jacobian(sys3, states), [states inputs], sys3_eq);
sys3_Blin = subs(jacobian(sys3, inputs), [states inputs], sys3_eq);

A3 = double(vpa(sys3_Alin));
B3 = double(vpa(sys3_Blin));


% Transfer Functions
syms s
sI = eye(12)*s;

% state space matrices
A = [0 0 0 1 zeros([8 1]).';
    0 0 0 0 1 zeros([7 1]).';
    0 0 0 0 0 1 zeros([6 1]).';
    0 0 0 -kx/m 0 0 g*sin(psi_d) g*cos(psi_d) 0 0 0 0;
    0 0 0 0 -ky/m 0 -g*cos(psi_d) g*sin(psi_d) 0 0 0 0;
    0 0 0 0 0 -kz/m 0 0 0 0 0 0;
    zeros([1 9]) 1 0 0;
    zeros([1 10]) 1 0;
    zeros([1 11]) 1;
    zeros([1 9]) -kp/Ix 0 0;
    zeros([1 10]) -kq/Iy 0;
    zeros([1 11]) -kr/Iz;];

B = [zeros([5 4]);
    1/m 0 0 0;
    zeros([3 4]);
    0 1/Ix 0 0;
    0 0 1/Iy 0;
    0 0 0 1/Iz];

C = [eye(3) zeros([3 9]);
    zeros([1 8]) 1 0 0 0];

D = zeros([4 4]);

% compute transfer function of system
inv = inv(sI - A);
Gs = C*inv*B + D;

% Analyse output of system
syms u1 u2 u3 u4;
Us = [u1 u2 u3 u4].';
Y = Gs*Us;
Yz = Y(3);
Yp = Y(4);









