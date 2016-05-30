function [ u ] = MPC_controller( sys, Q, R, S, N, umin, umax, uslope_min, uslope_max,  xmin, xmax, x, uprec, add_on_fail )
%MPC_CONTROLLER return the input as computed by an mpc aontroller
%   Detailed explanation goes here

options = optimset('quadprog');
options = optimset(options, 'LargeScale', 'off', 'Display' , 'off');

Ac = [];
for i = 1:N
    Ac = [Ac;sys.a^i];
end

Bc = [];
for i= 1:N
    line = [];
    for j = 1:N
        if(j>i)
            el = zeros(size(sys.b));
        else
            el = sys.a^(i-j) * sys.b;
        end
        line = [line el];
    end
    Bc = [Bc;line];
end

Qc = kron(eye(N-1), Q);
Qc = blkdiag(Qc, S);

Rc = kron(eye(N), R);

H = 2 * Bc' * Qc * Bc + Rc;
f = 2 * x' * Ac' * Qc * Bc;
umin_constr = umin * ones(N, 1);
umax_constr = umax * ones(N, 1);

flag = -1;
while flag < 0
    [A_state, b_state] = build_state_contraint_matrix(Ac, Bc, N, xmin, xmax, x);
    [A_slope, b_slope] = build_slope_contraint_matrix(sys.Ts, N, uslope_min, uslope_max, uprec);
    A_leq = [A_state; A_slope];
    b_leq = [b_state; b_slope];
    [A_leq, b_leq] = remove_infinite_constraints(A_leq, b_leq);
    
    [u, ~, flag] = quadprog(H, f, A_leq, b_leq, [], [], umin_constr, umax_constr, [], options);
    [xmin, xmax] = handle_error(flag, xmin, xmax, add_on_fail);
end
u = u(1);
end

function [xmin, xmax] = handle_error(flag, xminprec, xmaxprec, add_on_fail)
xmin = xminprec;
xmax = xmaxprec;
if(flag == -3)
    warning('Unbounded optimisation problem. Execution aborted.');
    xmax = xmaxprec + add_on_fail;
    xmin = xminprec - add_on_fail;
end
if(flag == -2)
    warning('Infeasible optimisation problem. Execution aborted.');
    xmax = xmaxprec + add_on_fail;
    xmin = xminprec - add_on_fail;
end
end

function [A, b] = build_state_contraint_matrix(Ac, Bc, N, x_min, x_max, x)
x_max_vect = repmat(x_max, N, 1);
A = Bc;
b = x_max_vect - Ac * x;

x_min_vect = repmat(x_min, N, 1);
A = [A; -Bc];
b = [b; - x_min_vect + Ac * x];
end

function [A, b] = build_slope_contraint_matrix(Ts, N, uslope_min, uslope_max, uprec)
A_max = eye(N) + tril(-ones(N), -1) + tril(ones(N), -2);
b_max = ones(N, 1) .* Ts .* uslope_max;
b_max(1) = b_max(1) + uprec;

b_min = - ones(N, 1) .* Ts .* uslope_min;
b_min(1) = b_min(1) - uprec;

A = [A_max; - A_max];
b = [b_max; b_min];
end

function [A, b] = remove_infinite_constraints(A, b)
A = A(b > -inf & b < inf, :);
b = b(b > -inf & b < inf);
end
