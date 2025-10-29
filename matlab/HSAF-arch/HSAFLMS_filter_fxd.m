function [error, s_vec, U_vec, U_vec_C, j_vec] = HSAFLMS_filter_fxd(x, s_vec, U_vec, U_vec_C, j_vec, d, w, q, par, fraclen, wrdlen)
    
    P = par.P;      
    Q = par.Q;
    delX = par.delX ;
    C = par.C;

    x_fi = x .* 2^(fraclen);
    x_fi = nearest(x_fi);

    d_fi = d .* 2^(fraclen);
    d_fi = nearest(d_fi);

    w_fi = w .* 2^(fraclen);
    w_fi = nearest(w_fi);

    q_fi = q .* 2^(fraclen);
    q_fi = nearest(q_fi);
    
%     for i=1:30
%     tmpvar(i) = sprintf("q_weight[%d] <= 16'h%s;",i-1,q_prnt(i,:)');
%     end

    s_vec_fi = s_vec .* 2^(fraclen);
    s_vec_fi = nearest(s_vec_fi);

    u_vec = zeros(P+1,1);
    u_vec_fi = zeros(P+1,1);
    u_vec(P+1) = 1;
    u_vec_fi(P+1) = 1*2^fraclen;

    x_div_delX_fi = sign(x_fi)*bitshift(abs(x_fi),log2(1/delX));
    if(overflow_det(x_div_delX_fi*2^(-fraclen),fraclen,wrdlen))
        x_div_delX_fi
    end

    if sign(x_div_delX_fi) == -1
        ceil_x_div_delX = sign(x_div_delX_fi)*bitshift(abs(x_div_delX_fi),-fraclen)-1;
    else
        ceil_x_div_delX = bitshift(x_div_delX_fi,-fraclen);
    end

    u_fi = x_div_delX_fi - ceil_x_div_delX*2^fraclen;
    u = u_fi * 2^(-fraclen);

    j = ceil_x_div_delX + (Q-1)/2;
        
    j_vec = [j, j_vec(1:end-1)];
    
    %u_vec = [u,1]';
    u_vec(P) = u;
    u_vec_fi(P) = u_fi;
    for p_ind = 2:P
%             u_vec = [u^p_ind; u_vec];
        u_vec_mult_fi = u_vec_fi(P-p_ind+2) * u_fi;
        u_vec_mult_fi_fraclen = fraclen + fraclen;
        u_vec_mult = u_vec_mult_fi * 2^(-u_vec_mult_fi_fraclen);
        u_vec_mult_fi_fraclen = fraclen;
        u_vec_mult_fi = u_vec_mult * 2^(u_vec_mult_fi_fraclen);
        u_vec_mult_fi = nearest(u_vec_mult_fi);
        u_vec_mult = u_vec_mult_fi * 2^(-u_vec_mult_fi_fraclen);
        u_vec_fi(P-p_ind+1) = u_vec_mult_fi;
        u_vec(P-p_ind+1) = u_vec_mult;
        
    end
%     end

    if(overflow_det(u_vec,fraclen,wrdlen))
        u_vec
    end

    U_vec = [u_vec, U_vec(:,1:end-1)];

    u_vec_C_fi = floor(u_vec_fi'*C);
    u_vec_C = u_vec_C_fi * 2^(-fraclen);
    U_vec_C = [u_vec_C', U_vec_C(:,1:end-1)];
    if(overflow_det(U_vec_C,fraclen,wrdlen))
        U_vec_C
    end

    s_mult_fi = u_vec_C_fi.*q_fi(j+1:j+P+1)';
    s_mult_fi_fraclen = fraclen + fraclen;
    s_mult = s_mult_fi * 2^(-s_mult_fi_fraclen);
    s_mult_fi_fraclen = fraclen;
    s_mult_fi = s_mult * 2^(s_mult_fi_fraclen);
    s_mult_fi = nearest(s_mult_fi);
    s_mult = s_mult_fi * 2^(-s_mult_fi_fraclen);
    if(overflow_det(s_mult,fraclen,wrdlen))
        s_mult
    end

    s_fi = sum(s_mult_fi);
    s_fraclen = fraclen;
    s = s_fi * 2^(-s_fraclen);
    s_fraclen = fraclen;
    s_fi = s * 2^(s_fraclen);
    s_fi = nearest(s_fi);
    if(overflow_det(s,fraclen,wrdlen))
        s
    end

    s_vec_fi = [s_fi, s_vec_fi(1:end-1)];
    s_vec = [s, s_vec(1:end-1)];

    y_mult_fi = w_fi .* [0*2^fraclen, s_vec_fi(1:length(w)-1)];
    y_mult_fraclen = fraclen + fraclen;
    y_mult = y_mult_fi * 2^(-y_mult_fraclen);
    y_mult_fraclen = fraclen;
    y_mult_fi = y_mult * 2^(y_mult_fraclen);
    y_mult_fi = nearest(y_mult_fi);
    if(overflow_det(y_mult,fraclen,wrdlen))
        y_mult
    end

    y_fi = sum(y_mult_fi);
    y_fraclen = fraclen;
    y = y_fi * 2^(-y_fraclen);
    y_fraclen = fraclen;
    y_fi = y * 2^(y_fraclen);
    y_fi = nearest(y_fi);
    if(overflow_det(y,fraclen,wrdlen))
        y
    end

    error_fi = d_fi - y_fi;
    error_fraclen = fraclen;
    error = error_fi * 2^(-error_fraclen);
    error_fraclen = fraclen;
    error_fi = error * 2^(error_fraclen);
    error_fi = nearest(error_fi);
    if(overflow_det(error,fraclen,wrdlen))
        error
    end
