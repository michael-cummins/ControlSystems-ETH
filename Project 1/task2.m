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

% Eigenvalues 
eigA = eig(A);
J = jordan(A)

% Answer = no
% Not all jordan blocks associate with eig = 0 are 1x1 
% Apparent from the [0 1; 0 0] block in the first entry











