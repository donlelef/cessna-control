% Stabilizing controller for Cessna Citation 500 aircraft
clear all;
close all;
warning('off', 'all');
clc;

%% Model definition
% x1: angle of attack
% x2: pitch angle
% x3: pitch rate
% x4: altitude
A = [-1.2822, 0, 0.98, 0; 0, 0, 1, 0; -5.4293, 0, -1.8366, 0; -128.2, 128.2, 0, 0];
B = [-0.3; 0; -17; 0];
C = [0, 1, 0, 0; 0, 0, 0, 1; -128.2, 128.2, 0, 0];
D = [0; 0; 0];
Ts = 0.1;
sys = ss(A, B, C, D);
sys_discrete = c2d(sys, Ts);

set_point = [0; 0; 0; 0];
T_sim = 10;
x0 = [0; 0; 0; 10];
N = 10;

%% Parameters
u_min = -0.262;
u_max= 0.262;
u_slew_rate_min = -0.524;
u_slew_rate_max = 0.524;
x2_min = -0.349;
x2_max = 0.349;
overshot_fraction_h = 1e-2;
overshot_fraction_l = 1e-2;
add_on_fail = overshot_fraction_h * [0, 0, 0, 10]';

%% 1. LQ controller
Q = eye(4);
R = 1;
[Kdlqr, S, closeLoopEigs] = dlqr(sys_discrete.a, sys_discrete.b, Q, R);
closeLoopMatrix = sys_discrete.a - sys_discrete.b * Kdlqr;
display(['Closed loop matrix A-BK eignevalues = ' num2str(closeLoopEigs')]);
if(abs(closeLoopEigs) < 1)
    display('Eigenvalues inside unit circle -> closed loop system is AS.')
else
    display('Eigenvalues outside unit circle -> error.')
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
display('Starting simulation without constraints.')
Q_simulink = Q
S_simulink = S
u_min_simulink = -inf
u_max_simulink = inf
x_min_simulink =  [-inf, -inf, -inf, -inf]'
x_max_simulink =  [+inf, +inf, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4. MPC controller with constraint on input variable
close all;
display('Starting simulation with constraints on input variable.')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, -inf, -inf, -inf]'
x_max_simulink =  [+inf, +inf, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4a. MPC more aggressive (Q = diag([1, 1, 1, 1e5]), S recomputed)
close all;
display('Starting simulation with more aggressive control on x4: Q = diag([1, 1, 1, 1e5]), S recomputed.')
Q_aggressive = diag([1, 1, 1, 1e5]);
[Kdlqr_aggressive, S_aggressive, ~] = dlqr(sys_discrete.a, sys_discrete.b, Q_aggressive, R);
Q_simulink = Q_aggressive
S_simulink = S_aggressive
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, -inf, -inf, -inf]'
x_max_simulink =  [+inf, +inf, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr_aggressive

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 4b. MPC less aggressive (Q = 1e-5 I, S recomputed)
close all;
display('Starting simulation where input is heavily penalised: Q = 1e5 * I, S recomputed.')
Q_aggressive = 1e-5 .* eye(4);
[Kdlqr_aggressive, S_aggressive, ~] = dlqr(sys_discrete.a, sys_discrete.b, Q_aggressive, R);
Q_simulink = Q_aggressive
S_simulink = S_aggressive
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, -inf, -inf, -inf]'
x_max_simulink =  [+inf, +inf, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr_aggressive

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 5. Addition of pitch angle constraint (state x2)
close all;
display('Starting simulation with constraints on input variable and pitch angle (x2).')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, -inf]'
x_max_simulink =  [+inf, x2_max, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 5b Removing input contraint to let x2 saturate.
close all;
display('Starting simulation with constraints on input variable and pitch angle (x2), no input constraint.')
Q_simulink = Q
S_simulink = S
u_min_simulink = -inf
u_max_simulink = inf
x_min_simulink =  [-inf, x2_min, -inf, -inf]'
x_max_simulink =  [+inf, x2_max, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 6. Addition of altitude overshot constraint (state x4)
close all;
display('Starting simulation with constraints on input variable, pitch angle (x2) and altitude overshot (x4).')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, - x0(4) * overshot_fraction_l]'
x_max_simulink =  [+inf, x2_max, +inf, x0(4) * (1 + overshot_fraction_h)]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;

%% 7. Addition of input slope constraint
display('Starting simulation with constraints on input variable, pitch angle (x2), altitude overshot (x4) and input slew rate.')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, -inf]'
x_max_simulink =  [+inf, x2_max, +inf, inf]'
uslopemin_simulink = u_slew_rate_min
uslopemax_simulink = u_slew_rate_max
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;
