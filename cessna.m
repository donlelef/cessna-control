% Stabilizing controller for Cessna Citation 500 aircraft
clear all;
close all;
clc;

%% Model definition
A = [-1.2822, 0, 0.98, 0; 0, 0, 1, 0; -5.4293, 0, -1.8366, 0; -128.2, 128.2, 0, 0];
B = [-0.3; 0; -17; 0];
C = [0, 1, 0, 0; 0, 0, 0, 1; -128.2, 128.2, 0, 0];
D = [0; 0; 0];
Ts = 0.1;
sys = ss(A, B, C, D);
sys_discrete = c2d(sys, Ts);

%% Parameters
u_min = -0.262;
u_max= 0.262;
u_slew_rate_min = -0.524;
u_slew_rate_max = 0.524;
x2_min = -0.349; % Pitch angle
x2_max = 0.349;

%% 1. LQ controller
Q = eye(4);
R = 1;
[Kdlqr, S, closeLoopEigs] = dlqr(sys_discrete.a, sys_discrete.b, Q, R);
closeLoopMatrix = sys_discrete.a - sys_discrete.b * Kdlqr;
display(['Closed loop matrix A-BK eignevalues = ' num2str(closeLoopEigs')]);
if(abs(closeLoopEigs) < 1)
    display('Eigenvalues inside unit circle -> closed loop system is AS.')
end

%% 2. Controller tuning with simulink model
set_point = [0; 0; 0; 0];
x0 = [0; 0; 0; 10];
u_min_simulink = -inf;
u_max_simulink = inf;
T_sim = 5;
sim('LQR_discrete');
open('LQR_discrete');
%pause;

u_min_simulink = u_min;
u_max_simulink = u_max;
T_sim = 10;
sim('LQR_discrete');
open('LQR_discrete');
%pause;

%% 3. MPC controller without active constraints
close all;
set_point = [0; 0; 0; 0];
N = 100;
x0 = [0; 0; 0; 10];
u_min_simulink = -inf;
u_max_simulink = inf;
u_mpc = MPC_controller(sys_discrete.a, sys_discrete.b, Q, R, S, N, u_min_simulink, u_max_simulink, x0);
T_sim = 10;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

