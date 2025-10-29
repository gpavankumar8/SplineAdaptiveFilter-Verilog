function [w_upd, q_upd] = HSAFLMS_update_fxd(s_vec, U_vec, U_vec_C, j_vec, error, w, q, par, fraclen, wrdlen)

    mu_w = par.mu_w_HSAFLMS;
    mu_q = par.mu_q_HSAFLMS;
    P =    par.P;
    j = j_vec(1);

    w_fi = w .* 2^(fraclen);
    w_fi = nearest(w_fi);

    q_fi = q .* 2^(fraclen);
    q_fi = nearest(q_fi);
    q_upd_fi = q_fi;

    error_fi = error .* 2^(fraclen);
    error_fi = nearest(error_fi);

    s_vec_fi = s_vec .* 2^(fraclen);
    s_vec_fi = nearest(s_vec_fi);

    U_vec_fi = U_vec .* 2^(fraclen);
    U_vec_fi = nearest(U_vec_fi);

    U_vec_C_fi = U_vec_C .* 2^(fraclen);
    U_vec_C_fi = nearest(U_vec_C_fi);


%     q(j+1:j+4) = q(j+1:j+4) + mu_q * error*C'*U_vec(:,j_vec==j)*w(j_vec==j)';

    error_mu_q_fi = 2^(nearest(log2(mu_q))) * error_fi;
    error_mu_q_fraclen = fraclen;
    error_mu_q = error_mu_q_fi * 2^(-error_mu_q_fraclen);
    error_mu_q_fraclen = fraclen;
    error_mu_q_fi = error_mu_q * 2^(error_mu_q_fraclen);
    error_mu_q_fi = nearest(error_mu_q_fi);
    if(overflow_det(error_mu_q,fraclen,wrdlen))
        error_mu_q
    end

    w_vmm_fi = w_fi(2:end);

    U_vec_C_w_mult_fi = U_vec_C_fi(:,j_vec==j) .* w_vmm_fi(j_vec==j);
    U_vec_C_w_mult_fraclen = fraclen + fraclen;
    U_vec_C_w_mult = U_vec_C_w_mult_fi * 2^(-U_vec_C_w_mult_fraclen);
    U_vec_C_w_mult_fraclen = fraclen;
    U_vec_C_w_mult_fi = U_vec_C_w_mult * 2^(U_vec_C_w_mult_fraclen);
    U_vec_C_w_mult_fi = nearest(U_vec_C_w_mult_fi);
    if(overflow_det(U_vec_C_w_mult,fraclen,wrdlen))
        U_vec_C_w_mult
    end

    U_vec_C_w_fi = sum(U_vec_C_w_mult_fi,2);
    U_vec_C_w_fraclen = fraclen;
    U_vec_C_w = U_vec_C_w_fi * 2^(-U_vec_C_w_fraclen);
    U_vec_C_w_fraclen = fraclen;
    U_vec_C_w_fi = U_vec_C_w * 2^(U_vec_C_w_fraclen);
    U_vec_C_w_fi = nearest(U_vec_C_w_fi);
    if(overflow_det(U_vec_C_w,fraclen,wrdlen))
        U_vec_C_w
    end

    error_mu_q_U_vec_C_w_fi = error_mu_q_fi .* U_vec_C_w_fi; 
    error_mu_q_U_vec_C_w_fi_fraclen = fraclen + fraclen;
    error_mu_q_U_vec_C_w = error_mu_q_U_vec_C_w_fi * 2^(-error_mu_q_U_vec_C_w_fi_fraclen);
    error_mu_q_U_vec_C_w_fi_fraclen = fraclen;
    error_mu_q_U_vec_C_w_fi = error_mu_q_U_vec_C_w * 2^(error_mu_q_U_vec_C_w_fi_fraclen);
    error_mu_q_U_vec_C_w_fi = nearest(error_mu_q_U_vec_C_w_fi);
    if(overflow_det(error_mu_q_U_vec_C_w,fraclen,wrdlen))
        error_mu_q_U_vec_C_w
    end

    q_upd_fi(j+1:j+P+1) = q_fi(j+1:j+P+1) + error_mu_q_U_vec_C_w_fi;
    q_upd_fraclen = fraclen;
    q_upd = q_upd_fi * 2^(-q_upd_fraclen);
    q_upd_fraclen = fraclen;
    q_upd_fi = q_upd * 2^(q_upd_fraclen);
    q_upd_fi = nearest(q_upd_fi);
    q_upd = q_upd_fi * 2^(-q_upd_fraclen);
    if(overflow_det(q_upd,fraclen,wrdlen))
        q_upd
    end

    % W Update

    error_mu_w_fi = 2^(nearest(log2(mu_w))) * error_fi;
    error_mu_w_fraclen = fraclen;
    error_mu_w = error_mu_w_fi * 2^(-error_mu_w_fraclen);
    error_mu_w_fraclen = fraclen;
    error_mu_w_fi = error_mu_w * 2^(error_mu_w_fraclen);
    error_mu_w_fi = nearest(error_mu_w_fi);
    if(overflow_det(error_mu_w,fraclen,wrdlen))
        error_mu_w
    end

    s_vec_mu_w_err_fi = [0*2^fraclen, s_vec_fi(1:length(w)-1)] * error_mu_w_fi;
    s_vec_mu_w_err_fraclen = fraclen + error_mu_w_fraclen;
    s_vec_mu_w_err = s_vec_mu_w_err_fi .* 2^(-s_vec_mu_w_err_fraclen);
    s_vec_mu_w_err_fraclen = fraclen;
    s_vec_mu_w_err_fi = s_vec_mu_w_err .* 2^(s_vec_mu_w_err_fraclen);
    s_vec_mu_w_err_fi = nearest(s_vec_mu_w_err_fi);
    if(overflow_det(s_vec_mu_w_err,fraclen,wrdlen))
        s_vec_mu_w_err
    end

    w_upd_fi = w_fi + s_vec_mu_w_err_fi;
    w_upd_fraclen = fraclen;
    w_upd = w_upd_fi * 2^(-w_upd_fraclen);
    w_upd_fraclen = fraclen;
    w_upd_fi = w_upd * 2^(w_upd_fraclen);
    w_upd_fi = nearest(w_upd_fi);
    w_upd = w_upd_fi * 2^(-w_upd_fraclen);
    if(overflow_det(w_upd,fraclen,wrdlen))
        w_upd
    end