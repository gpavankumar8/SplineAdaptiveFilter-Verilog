function [error, w] = RFFKLMS_Gaussian_filter(x_vec, d, w, par)

    sigma = par.sigma;
    D = par.D;
    mu = par.muRFFKLMS;

    R = par.R;
    b = par.b;

%     z= sqrt(2/D)*cos(R*x_vec' + b);
    z = exp(-(vecnorm((repmat(x_vec,[length(R),1])-R),2,2).^2)./(2*b.^2));
    error = d - w' * [0;z]; 
    w = w + mu * error * [0;z];