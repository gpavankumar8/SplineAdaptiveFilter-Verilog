function [error, y] = LMLS_filter(x_vec, d, w)

    y = [0,x_vec] * w';
%     y = x_vec * w';
    error = d - y;
    
end