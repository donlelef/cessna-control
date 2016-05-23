% Stabilizing controller for Cessna Citation 500 aircraft
clear all;
close all;
warning('off', 'all');
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
x_min_simulink =  [-inf, -inf, -inf, -inf]';
x_max_simulink =  [+inf, +inf, +inf, +inf]';
u_min = -0.262;
u_max= 0.262;
u_slew_rate_min = -0.524;
u_slew_rate_max = 0.524;
x2_min = -0.349; % Pitch angle
x2_max = 0.349;
overshot_fraction_h = 1e-10;
overshot_fraction_l = 1e-10;
set_point = [0; 0; 0; 0];
T_sim = 10;
x0 = [0; 0; 0; 10];
N = 10;

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
u_min_simulink = -inf;
u_max_simulink = inf;
Kdlqr_simulink = Kdlqr;
sim('LQR_discrete');
open('LQR_discrete');
%pause;

%% 3. MPC controller without active constraints
close all;
Q_simulink = Q;
S_simulink = S;
u_min_simulink = -inf;
u_max_simulink = inf;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4. MPC controller with constraint on input variable
close all;
u_min_simulink = u_min;
u_max_simulink = u_max;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4a. MPC more aggressive (Q = 100I, S recomputed)
close all;
Q_aggressive = 100 .* Q;
[Kdlqr_aggressive, S_aggressive, ~] = dlqr(sys_discrete.a, sys_discrete.b, Q_aggressive, R);
Q_simulink = Q_aggressive;
S_simulink = S_aggressive;
Kdlqr_simulink = Kdlqr_aggressive;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4b. MPC more aggressive (Q = 100 I, S computed with Q = I)
close all;
S_simulink = S;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4c. MPC more aggressive (Q = I, S computed with Q = 100 I)
close all;
Q_simulink = Q;
S_simulink = S_aggressive;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 5. Addition of pitch angle constraint (state x2)
close all;
S_simulink = S;
x_min_simulink(2) = x2_min;
x_max_simulink(2) = x2_max;
Kdlqr_simulink = Kdlqr;
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 6. Addition of altitude overshot constraint (state x4)
close all;
x_min_simulink(4) = - x0(4) * overshot_fraction_h;
x_max_simulink(4) = x0(4) * (1 + overshot_fraction_l);
sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;