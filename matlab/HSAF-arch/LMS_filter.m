function [error, y] = LMS_filter(x_vec, d, w)

    y = [0,x_vec] * w';
%     y = x_vec * w';
    error = d - y;
