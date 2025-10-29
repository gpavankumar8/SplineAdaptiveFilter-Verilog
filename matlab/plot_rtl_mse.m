
monte = 1;
N = 20000;

fpath = '../DHSAF-I/simulation/outputs/L8_v5/';
% fpath = '../DHSAF-I-CG/simulation/outputs/L8_v5/';
% fpath = '../DHSAF-II/simulation/outputs/L8_v5/';
% fpath = '../DHSAF-II-CG/simulation/outputs/L8_v5/';
% fpath = '../HSAF/simulation/outputs/L8_v5/';


fraclen = 12;
wrdlen = 16;

error_mc_rtl = zeros(N,1);

for trial=1:monte
    
    err_rtl = rtl_bin2dec_err(sprintf([fpath 'error_rtl%i.txt'],trial), fraclen, wrdlen);
    error_mc_rtl = error_mc_rtl + (err_rtl.^2);

end

error_mc_rtl = movmean(error_mc_rtl / monte,500);

figure;
plot(10*log10(error_mc_rtl));
grid on
xlabel('Iteration');
ylabel('MSE (dB)');
title('Adaptive filter convergence');

steadystate_MSE = mean(10*log10(error_mc_rtl(end-1000:end)));
disp("Steady state MSE");
disp(steadystate_MSE);