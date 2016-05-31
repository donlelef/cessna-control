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

%% Open Loop Analysis
pzmap(sys_discrete);

%% 1. LQ controller
display('1. LQ controller synthesis');
close all;
Q = eye(4)
R = 1
[Kdlqr, S, closeLoopEigs] = dlqr(sys_discrete.a, sys_discrete.b, Q, R);
theorem = (rank(ctrb(sys_discrete)) == 4) && (all(eig(S)>0));
if(theorem)
    display('The assumptions for the asymptotical stability of LQR are verified');
end
display(['Closed loop matrix A-BK eignevalues = ' num2str(closeLoopEigs')]);
if(abs(closeLoopEigs) < 1)
    display('Eigenvalues inside unit circle -> closed loop system is AS.')
else
    display('Eigenvalues outside unit circle -> error.')
end

%% 2. Controller tuning with simulink model
display('2. LQ controller simulation');
u_min_simulink = -inf
u_max_simulink = inf
Kdlqr_simulink = Kdlqr
sim('LQR_discrete');
open('LQR_discrete');
%pause;

%% 3. MPC controller without active constraints
display('3. MPC addition and simulation without constraints');
close all;
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
display('4. Simulation with constraint on u');
close all;
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
display('4a. Simulation with more aggressive control on x4: Q = diag([1, 1, 1, 1e5]), S recomputed.');
close all;
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
display('Simulation where input is heavily penalised: Q = 1e-5 * I, S recomputed.');
close all;
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
display('Addition of constraint on x2 within the MPC optimization')
close all;
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

%% 5a Removing input contraint to let x2 saturate.
display('5a. Simulation with no constraint on u and constraints on x2')
close all;
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
display('6. Simulation with constraints on u, x2, and x4 overshot');
close all;
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
display('Simulation with constraints on u, x2, x4 overshot and u slew rate')
close all;
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, - x0(4) * 1 *overshot_fraction_l]'
x_max_simulink =  [+inf, x2_max, +inf, x0(4) * (1 + overshot_fraction_h)]'
uslopemin_simulink = u_slew_rate_min
uslopemax_simulink = u_slew_rate_max
Kdlqr_simulink = Kdlqr

sim('MPC_vs_LQR');
open('MPC_vs_LQR');
% pause;
