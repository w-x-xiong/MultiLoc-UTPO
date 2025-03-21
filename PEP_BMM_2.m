function [y_breve_tilde, iter_num] = PEP_BMM_2(Rx, Rg, dRg, SigmaRg, SigmadRg, y_ini, thre, NM)

scale = 1;
Rx = Rx/scale;
Rg = Rg/scale;
dRg = dRg/scale;

[~, L] = size(Rg);
[H, ~] = size(Rx);

iter_num = 0;

y_breve_tilde = y_ini/scale;

while iter_num <= (NM-1)
    y_breve_tilde_old = y_breve_tilde;
    mu = zeros(1, L);
    xi = zeros(1, L);
    phi = zeros(1, L);
    zeta = zeros(H, L);
    eta = zeros(H, L);
    rho = zeros(1, L);
    tau = zeros(1, L);
    for l = 1:L
        mu(1,l) = Rg(1,l) - norm(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H)) - norm(y_breve_tilde(1:H) - Rx(:,l)) + y_breve_tilde(2*H+1);
        xi(1,l) = mu(1,l)/(2*SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H)));
        phi(1,l) = mu(1,l)/(2*SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - Rx(:,l)));
        zeta(:,l) = -2*Rg(1,l)*(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H))/(SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H)));
        eta(:,l) = -2*Rg(1,l)*(y_breve_tilde(1:H) - Rx(:,l))/(SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - Rx(:,l)));
        rho(1,l) = norm(y_breve_tilde(1:H) - Rx(:,l))/(SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H)));
        tau(1,l) = norm(y_breve_tilde(1:H) - y_breve_tilde(H+1:2*H))/(SigmaRg(1, l)^2*norm(y_breve_tilde(1:H) - Rx(:,l)));
    end
    x_subp_numerator = zeros(H,1);
    x_subp_denominator = 0;
    for l = 1:L
        x_subp_numerator = x_subp_numerator + 2*(xi(1,l) + 1/SigmaRg(1, l)^2 + rho(1,l))*y_breve_tilde(H+1:2*H) + 2*(1/SigmaRg(1, l)^2 + tau(1,l) + phi(1,l))*Rx(:,l) - zeta(:,l) - eta(:,l);
        x_subp_denominator = x_subp_denominator + 2*(xi(1,l) + 2/SigmaRg(1, l)^2 + rho(1,l) + tau(1,l) + phi(1,l));
    end
    y_breve_tilde(1:H) = x_subp_numerator/x_subp_denominator;

    omega0_subp_numerator = 0;
    omega0_subp_denominator = 0;
    for l = 1:L
        omega0_subp_numerator = omega0_subp_numerator + mu(1,l)/SigmaRg(1, l)^2 - 2*(norm(y_breve_tilde(H+1:2*H) - Rx(:,l)) - dRg(1,l))/SigmadRg(1, l)^2;
        omega0_subp_denominator = omega0_subp_denominator + 2*(1/SigmaRg(1, l)^2 + 1/SigmadRg(1, l)^2);
    end
    y_breve_tilde(2*H+1) = omega0_subp_numerator/omega0_subp_denominator;

    nu = zeros(H, L);
    psi = zeros(1, L);
    epsilon = zeros(1, L);
    kappa = zeros(H, L);
    for l = 1:L
        nu(:,l) = -2*Rg(1,l)*(y_breve_tilde(H+1:2*H) - y_breve_tilde(1:H))/(SigmaRg(1, l)^2*norm(y_breve_tilde(H+1:2*H) - y_breve_tilde(1:H)));
        psi(1,l) = (y_breve_tilde(2*H+1)+norm(y_breve_tilde(1:H) - Rx(:,l)))/(SigmaRg(1, l)^2*norm(y_breve_tilde(H+1:2*H) - y_breve_tilde(1:H)));
        epsilon(1,l) = y_breve_tilde(2*H+1)/(SigmadRg(1, l)^2*norm(y_breve_tilde(H+1:2*H) - Rx(:,l)));
        kappa(:,l) = -2*dRg(1,l)*(y_breve_tilde(H+1:2*H) - Rx(:,l))/(SigmadRg(1, l)^2*norm(y_breve_tilde(H+1:2*H) - Rx(:,l)));
    end
    t_subp_numerator = zeros(H,1);
    t_subp_denominator = 0;
    for l = 1:L
        t_subp_numerator = t_subp_numerator + 2*y_breve_tilde(1:H)*(1/SigmaRg(1, l)^2 + psi(1,l)) + 2*Rx(:,l)*(1/SigmadRg(1, l)^2 + epsilon(1,l)) - (nu(:,l) + kappa(:,l));
        t_subp_denominator = t_subp_denominator + 2/SigmaRg(1, l)^2 + 2*psi(1,l) + 2/SigmadRg(1, l)^2 + 2*epsilon(1,l);
    end
    y_breve_tilde(H+1:2*H) = t_subp_numerator/t_subp_denominator;
    
    iter_num = iter_num + 1;
    if (norm(y_breve_tilde(1:H) - y_breve_tilde_old(1:H))/norm(y_breve_tilde_old(1:H)) <= thre)
        break
    end
end

y_breve_tilde = y_breve_tilde*scale;

end

