function [u, Tx, theta_est] = WLS_XWu_2(Rx, Rg, dRg, Sigma_alpha)

[M, N] = size(Rg);
[p, ~] = size(Rx);

Gr = zeros(M*N, (p*(M+1)+4*M));
hr = zeros(M*N, 1);
Gd = zeros(M*N, (p*(M+1)+4*M));
hd = zeros(M*N, 1);

for i = 1:M
    for j = 1:N
        k = (i-1)*N+j;
        Gr(k, 1:p) = -Rx(:,j)';
        Gr(k, p*(M+1)+i) = Rg(i,j);
        Gr(k, p*(M+1)+M+3*i+(-2:0)) = [Rg(i,j), 1, 0.5];
        hr(k) = 0.5*(Rg(i,j)^2 - norm(Rx(:,j))^2);
        Gd(k, i*p+(1:p)) = -Rx(:,j)';
        Gd(k, p*(M+1)+i) = dRg(i,j);
        Gd(k, p*(M+1)+M+3*i) = 0.5;
        hd(k) = 0.5*(dRg(i,j)^2 - norm(Rx(:,j))^2);
    end
end

h1 = [hr', hd']';
G1 = [Gr', Gd']';

%Stage-one solution - coarse
Psi = (pinv(G1'*eye(2*M*N)*G1))*G1'*eye(2*M*N)*h1;

u = Psi(1:p);

Bri = zeros(N, N, M);

for i = 1:M
    for j = 1:N
        Bri(j, j, i) = norm(u - Rx(:,j));
    end
end

Br = [];

for i = 1:M
    Br = blkdiag(Br, Bri(:, :, i));
end

Tx = [];

for i = 1:M
    Tx = [Tx, Psi(i*p+(1:p))];
end

Bdi = zeros(N, N, M);

for i = 1:M
    for j = 1:N
        Bdi(j, j, i) = norm(Tx(:,i) - Rx(:,j));
    end
end

Bd = [];

for i = 1:M
    Bd = blkdiag(Bd, Bdi(:, :, i));
end

Pri = zeros(N, p*N, M);
Pdi = zeros(N, p*N, M);

for i = 1:M
    for j = 1:N
        Pri(j, (j-1)*p+1:j*p, i) = ((Rx(:,j) - u)/norm(Rx(:,j) - u))';
        Pdi(j, (j-1)*p+1:j*p, i) = ((Rx(:,j) - Tx(:,i))/norm(Rx(:,j) - Tx(:,i)))';
    end
end

Pr = [];
Pd = [];

for i = 1:M
    Pr = [Pr,Pri(:,:,i)'];
    Pd = [Pd,Pdi(:,:,i)'];
end

Pr = Pr';
Pd = Pd';

B1 = [Br*Pr, Br, zeros(M*N,M*N); Bd*Pd, zeros(M*N,M*N), Bd];

W1 = pinv((B1*Sigma_alpha*B1'));

%Stage-one solution - refined
Psi = (pinv(G1'*W1*G1))*G1'*W1*h1;

%Stage-two solution
K = p*(M+1)+M;
G2 = zeros(K+3*M, K);
h2 = zeros(K+3*M, 1);
B2 = zeros(K+3*M, K+3*M);

u = Psi(1:p);

Tx = [];
theta = zeros(M,1);

for i = 1:M
    Tx = [Tx, Psi(i*p+(1:p))];
    theta(i) = Psi(p*(M+1)+i);
end

G2(1:K,1:K) = eye(K);

B2(1:K,1:K) = eye(K);

for i = 1:M
    yi = Psi(p*(M+1)+M+((i-1)*p+1:i*p));
    G2(K+i,1:p) = (u - Tx(:,i))';
    G2(K+i,i*p+(1:p)) = -(u - Tx(:,i))';
    G2(K+M+i,1:p) = Tx(:,i)';
    G2(K+M+i,i*p+(1:p)) = (u - 2*Tx(:,i))';
    G2(K+M+i,p*(M+1)+i) = -2*yi(1);
    G2(K+2*M+i,i*p+(1:p)) = Tx(:,i)';
    G2(K+2*M+i,p*(M+1)+i) = -theta(i);
    h2(1:K) = Psi(1:K);
    h2(K+i) = Psi(K+3*i-2)^2;
    h2(K+M+i) = 2*Psi(K+3*i-1);
    h2(K+2*M+i) = Psi(K+3*i);
    B2(K+i,1:p) = -(u - Tx(:,i))';
    B2(K+i,i*p+(1:p)) = (u - Tx(:,i))';
    B2(K+i,K+3*i-2) = 2*yi(1);
    B2(K+M+i,1:p) = -Tx(:,i)';
    B2(K+M+i,i*p+(1:p)) = -(u - 2*Tx(:,i))';
    B2(K+M+i,K+3*i-2) = 2*theta(i);
    B2(K+M+i,K+3*i-1) = 2;
    B2(K+2*M+i,i*p+(1:p)) = -Tx(:,i)';
    B2(K+2*M+i,p*(M+1)+i) = theta(i);
    B2(K+2*M+i,K+3*i) = 1;
end

W2 = pinv(B2*(pinv(G1'*W1*G1))*B2');

Phi = (pinv(G2'*W2*G2))*G2'*W2*h2;

u = Phi(1:p);

Tx = [];

for i = 1:M
    Tx = [Tx, Phi(i*p+(1:p))];
end

theta_est = [];

for i = 1:M
    theta_est = [theta_est; Phi(p*(M+1)+i)];
end

end

