function [w, q] = HSAFLMS_update(s_vec, U_vec, j_vec, error, w, q, par)

    mu_w = par.mu_w_HSAFLMS;
    mu_q = par.mu_q_HSAFLMS;
    C = par.C;
    P = par.P;

    j = j_vec(1);
%     j_vec == j
    w_vmm = w(2:end);

    q(j+1:j+P+1) = q(j+1:j+P+1) + mu_q * error*C'*U_vec(:,j_vec==j)*w_vmm(j_vec==j)';
    w = w + mu_w * error*[0,s_vec];
%     w = w + mu_w * error*[1,s_vec];