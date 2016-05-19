function [ u ] = MPC_controller( A, B, Q, R, S, N, umin, umax, x )
%MPC_CONTROLLER return the input as computed by an mpc aontroller
%   Detailed explanation goes here

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

H = Bc' * Qc' * Bc + Rc;
f = 2 * x' * Ac' * Qc * Bc;
A_constr = zeros(1, N);
b_constr = 0;
umin_constr = umin * ones(N, 1);
umax_constr = umax * ones(N, 1);
u = quadprog(H, f, A_constr, b_constr, A_constr, b_constr, umin_constr, umax_constr);
u = u(1);
end