function [x] = SDP_UTPO(Rx, Rg, dRg, Q)

[k, N] = size(Rx);

br_vec = [];
bd_vec = [];
Ar_mtx = [];
Ad_mtx = [];

for j = 1:N
    br_vec = [br_vec;Rg(1,j)^2 - norm(Rx(:,j))^2];
    bd_vec = [bd_vec;dRg(1,j)^2 - norm(Rx(:,j))^2];
    ar_vec = [-Rx(:,j)',zeros(1,k),Rg(1,j),Rg(1,j),-0.5,0]';
    ad_vec = [zeros(1,k),-Rx(:,j)',dRg(1,j),zeros(1,2),-0.5]';
    Ar_mtx = [Ar_mtx;ar_vec'];
    Ad_mtx = [Ad_mtx;ad_vec'];
end

br_vec = br_vec/2;
bd_vec = bd_vec/2;

b_vec = [br_vec',bd_vec']';

A_mtx = [Ar_mtx',Ad_mtx']';

W_mtx = pinv(Q);

Psi_mtx = [A_mtx'*W_mtx*A_mtx, -A_mtx'*W_mtx*b_vec; -b_vec'*W_mtx*A_mtx, b_vec'*W_mtx*b_vec];

cvx_begin quiet

variables x(2*k+4) X(2*k+4, 2*k+4)

expression obj

obj = trace([X,x;x',1]*Psi_mtx);

minimize obj

subject to

x(2*k+3) == trace(X((k+1):2*k,(k+1):2*k)) - 2*trace(X(1:k,(k+1):2*k)) + X(2*k+1,2*k+1) + 2*X(2*k+1,2*k+2);

x(2*k+4) == X(2*k+1,2*k+1) - trace(X((k+1):2*k,(k+1):2*k));

X(2*k+2,2*k+2) == trace(X(1:k,1:k)) - 2*trace(X(1:k,(k+1):2*k)) + trace(X((k+1):2*k,(k+1):2*k));

[X,x;x',1] == semidefinite(2*k+5);

cvx_end

diag_Br_elem_vec = [];
diag_Bd_elem_vec = [];

for j = 1:N
    diag_Br_elem_vec = [diag_Br_elem_vec;norm(x(1:k) - Rx(:,j))];
    diag_Bd_elem_vec = [diag_Bd_elem_vec;norm(x((k+1):2*k) - Rx(:,j))];
end

Br_mtx = diag(diag_Br_elem_vec);
Bd_mtx = diag(diag_Bd_elem_vec);

B_mtx = blkdiag(Br_mtx,Bd_mtx);

W_mtx = pinv((B_mtx*Q*B_mtx'));

Psi_mtx = [A_mtx'*W_mtx*A_mtx, -A_mtx'*W_mtx*b_vec; -b_vec'*W_mtx*A_mtx, b_vec'*W_mtx*b_vec];

cvx_begin quiet

variables x(2*k+4) X(2*k+4, 2*k+4)

expression obj

obj = trace([X,x;x',1]*Psi_mtx);

minimize obj

subject to

x(2*k+3) == trace(X((k+1):2*k,(k+1):2*k)) - 2*trace(X(1:k,(k+1):2*k)) + X(2*k+1,2*k+1) + 2*X(2*k+1,2*k+2);

x(2*k+4) == X(2*k+1,2*k+1) - trace(X((k+1):2*k,(k+1):2*k));

X(2*k+2,2*k+2) == trace(X(1:k,1:k)) - 2*trace(X(1:k,(k+1):2*k)) + trace(X((k+1):2*k,(k+1):2*k));

[X,x;x',1] == semidefinite(2*k+5);

cvx_end

end

