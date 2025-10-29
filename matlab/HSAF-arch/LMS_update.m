function [w] = LMS_update(x_vec, w, mu, error)

    w = w + mu * [0,x_vec] * error;
%      w = w + mu * x_vec * error;