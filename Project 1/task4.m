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
k = [-0.0001 0.005 0.00005 0];
x0 = [-1; 0];
t = (0:0.01:300);

% Convergent system oscillations
Af = [0 1; -k(1)/m -kz/m];
Bf = [0;0];
C = [1 0];
D = 0;
sys_feedback = ss(Af, Bf, C, D);
ya = initial(sys_feedback, x0, t);

% divergent system
Af = [0 1; -k(2)/m -kz/m];
sys_feedback = ss(Af, Bf, C, D);
yb = initial(sys_feedback, x0, t);

% convergent with oscillations
Af = [0 1; -k(3)/m -kz/m];
sys_feedback = ss(Af, Bf, C, D);
yc = initial(sys_feedback, x0, t);

% Convergent to non-zero value
Af = [0 1; -k(4)/m -kz/m];
sys_feedback = ss(Af, Bf, C, D);
yd = initial(sys_feedback, x0, t);

% Plotting
subplot(2,2,1)
plot(t, ya)
hold on
plot(t, zeros([1, length(t)]), 'r--')
title('Divergent, K = -0.0001')
xlabel('Time')
ylabel('Output')
axis([0 60 -2.5 0.5])
hold off

subplot(2,2,2)
plot(t, yb)
hold on
plot(t, zeros([1, length(t)]), 'r--')
title('Convergent with oscillations, K = 0.005')
xlabel('Time')
ylabel('Output')
axis([0 60 -1.5 1])
hold off

subplot(2,2,4)
plot(t, yc)
hold on
plot(t, zeros([1, length(t)]), 'r--')
title('Convergent without Oscillations, K = 0.00005')
xlabel('Time')
ylabel('Output')
axis([0 300 -1.5 1.5])
hold off


subplot(2,2,3)
plot(t, yd)
hold on
plot(t, zeros([1, length(t)]), 'r--')
title('Convergent to non-zero state, K = 0')
xlabel('Time')
ylabel('Output')
axis([0 300 -1.5 1.5])
hold off



