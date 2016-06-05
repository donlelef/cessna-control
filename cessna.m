% Stabilizing controller for Cessna Citation 500 aircraft
clear all;
close all;
warning('off', 'all');
clc

%% Model definition
% x1: angle of attack
% x2: pitch angle
% x3: pitch rate
% x4: altitude
% y1: putch angle
% y2: altitude
% y3: vertical speed
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
overshot_fraction_h = 2e-2;
overshot_fraction_l = 2e-2;

%% Open Loop Analysis
if(rank(ctrb(sys_discrete)) == 4 && rank(ctrb(sys_discrete)) == 4)
    display('Rank of controllability and observability matrix is 4 -> the system is controllable and observable');
end
pzmap(sys_discrete);

%% State observer
if(rank(obsv(sys)) == 4)
    display('The assumptions for the a1symptotical stability of state observer are verified');
end
poles = [0.01, 0.012, 0.014, 0.016];
L = place(sys_discrete.a', sys_discrete.c', poles)';
sys_observer = ss(sys_discrete.a - L*sys_discrete.c, [sys_discrete.b, L], eye(4), zeros(4, 4), Ts);
sys_observer_discrete = sys_observer;
display(['Closed loop matrix A-LC eignevalues = ' num2str(eig(sys_observer)')]);
if(all(abs(eig(sys_observer)) < 1))
    display('Closed loop matrix A-LC is AS.');
end
pzmap(sys_observer);

%% 1. LQ controller
display('1. LQ controller synthesis');
close all;
Q = eye(4)
R = 10
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

%% 2b. Adding input saturation to LQR regulator
display('2. Adding input saturation to LQR regulator');
u_min_simulink = u_min
u_max_simulink = u_max
Kdlqr_simulink = Kdlqr
sim('LQR_discrete');
open('LQR_discrete');

%% 3. MPC controller without active constraints
display('3. MPC addition and simulation without constraints');
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

%% 4a. MPC more aggressive (Q = diag([1, 1, 1, 100]), S recomputed)
display('4a. Simulation with more aggressive control on x4: Q = diag([1, 1, 1, 100]), S recomputed.');
close all;
Q_aggressive = diag([1, 1, 1, 10]);
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

%% 4b. MPC less aggressive (Q = 1e-5 I, S recomputed)
display('Simulation where input is heavily penalised: Q = 1e-5 * I, S recomputed.');
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

%% 5. Addition of pitch angle constraint (state x2)
display('Addition of constraint on x2 within the MPC optimization')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, -inf]'
x_max_simulink =  [+inf, x2_max, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('Complete_schema');
open('Complete_schema');

%% 5a Removing input contraint to let x2 saturate.
display('5a. Simulation with no constraint on u and constraints on x2')
Q_simulink = Q
S_simulink = S
u_min_simulink = -inf
u_max_simulink = inf
x_min_simulink =  [-inf, x2_min, -inf, -inf]'
x_max_simulink =  [+inf, x2_max, +inf, +inf]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('Complete_schema');
open('Complete_schema');

%% 6. Addition of altitude overshot constraint (state x4)
display('6. Simulation with constraints on u, x2, and x4 overshot');
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, - x0(4) * overshot_fraction_l]'
x_max_simulink =  [+inf, x2_max, +inf, x0(4) * (1 + overshot_fraction_h)]'
uslopemin_simulink = -inf
uslopemax_simulink = inf
Kdlqr_simulink = Kdlqr

sim('Complete_schema');
open('Complete_schema');

%% 7. Addition of input slope constraint
display('Simulation with constraints on u, x2, x4 overshot and u slew rate')
Q_simulink = Q
S_simulink = S
u_min_simulink = u_min
u_max_simulink = u_max
x_min_simulink =  [-inf, x2_min, -inf, - inf]'
x_max_simulink =  [+inf, x2_max, +inf, inf]'
uslopemin_simulink = u_slew_rate_min
uslopemax_simulink = u_slew_rate_max
Kdlqr_simulink = Kdlqr

sim('Complete_schema');
open('Complete_schema');
