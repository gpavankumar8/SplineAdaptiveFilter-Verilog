function [error, w] = RFFKLMS_filter(x_vec, d, w, par)

    sigma = par.sigma;
    D = par.D;
    mu = par.muRFFKLMS;

    R = par.R;
    b = par.b;

    z= sqrt(2/D)*cos(R*x_vec' + b);
    error = d - w' * [0;z]; 
    w = w + mu * error * [0;z];