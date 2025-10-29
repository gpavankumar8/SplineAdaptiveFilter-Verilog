function [err] = rtl_bin2dec_err(filename, fraclen, wrdlen)

fid = fopen(filename,'r');

formatspec = '%f';
A = fscanf(fid,formatspec);
A(1) = [];
q = num2str(A);
m = quantizer([wrdlen,fraclen]);
err = bin2num(m,q);

fclose(fid);