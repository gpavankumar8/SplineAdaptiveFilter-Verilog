function [w] = LMLS_update(x_vec, w, mu, alpha, error)

    w = w + mu * [0,x_vec] * (alpha*error^3)/(1+alpha*error^2);
%      w = w + mu * x_vec * error;

end