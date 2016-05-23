function [ u ] = MPC_controller( A, B, Q, R, S, N, umin, umax, xmin, xmax, x )
%MPC_CONTROLLER return the input as computed by an mpc aontroller
%   Detailed explanation goes here

options = optimset('quadprog');
options = optimset(options, 'LargeScale', 'off', 'Display' , 'off');

Ac = [];
for i = 1:N
    Ac = [Ac;A^i];
end

Bc = [];
for i= 1:N
    line = [];
    for j = 1:N
        if(j>i)
            el = zeros(size(B));
        else
            el = A^(i-j) * B;
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

[A_leq, b_leq] = build_inequality(Ac, Bc, N, xmin, xmax, x);

u = quadprog(H, f, A_leq, b_leq, [], [], umin_constr, umax_constr, [], options);
u = u(1);
end

function [A, b] = build_inequality(Ac, Bc, N, x_min, x_max, x)
x_max_vect = repmat(x_max, N, 1);
A = Bc;
b = x_max_vect - Ac * x;

x_min_vect = repmat(x_min, N, 1);
A = [A; -Bc];
b = [b; - x_min_vect + Ac * x];

A = A(b > -inf & b < inf, :);
b = b(b > -inf & b < inf);

end
