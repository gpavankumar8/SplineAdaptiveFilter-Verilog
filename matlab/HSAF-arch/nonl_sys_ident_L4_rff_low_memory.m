%clear all;

tic
rng(1);

% Parameters - experiment
monte = 50;    % No. of experiments
N = 20000;      % Number of samples
L = 4;         % Filter order
nv = 0.01;     % Noise variance

% Fixed point
fraclen = 12;
wrdlen = 16;

% Delayed update
del = 0;

% Parameters - LMS
mu_LMS = 0.05;       % Step size (1/64)

% RFFKLMS  parameters
RFF_par.muRFFKLMS = 0.025;
RFF_par.D=512;              
RFF_par.sigma = 5;
RFF_par.R = mvnrnd(zeros(L,1),(1/(RFF_par.sigma^2))*eye(L),RFF_par.D);
RFF_par.b = 2*pi*rand(RFF_par.D,1);


% WSAFLMS parameters
WSAFLMS_par.mu_w_WSAFLMS = 0.0105;       % Step size
WSAFLMS_par.mu_q_WSAFLMS = 0.0140;       % Step size
    
WSAFLMS_par.P = 3;      
WSAFLMS_par.Q = 25;
WSAFLMS_par.ip_range = 6;
WSAFLMS_par.delX = WSAFLMS_par.ip_range/(WSAFLMS_par.Q-1);
WSAFLMS_par.C = ...
    0.5*[-1, 3,-3, 1; 
          2,-5, 4,-1;
         -1, 0, 1, 0;
          0, 2, 0, 0];

% HSAFLMS parameters
HSAFLMS_par.mu_w_HSAFLMS = 1/64;       % Step size
HSAFLMS_par.mu_q_HSAFLMS = 1/64;       % Step size
HSAFLMS_par.P = WSAFLMS_par.P; 
HSAFLMS_par.Q = WSAFLMS_par.Q;
HSAFLMS_par.ip_range = WSAFLMS_par.ip_range;
HSAFLMS_par.delX = WSAFLMS_par.delX ;
HSAFLMS_par.C = WSAFLMS_par.C;


% TFLAF parameters
TFLAF_par.mu = 1/1024;     % Step size = 1/128
P_t = 3;
TFLAF_par.Q_t = 2*P_t+1;
TFLAF_par.expL = L*TFLAF_par.Q_t;

% RCTFLAF parameters
RCTFLAF_par.mu_w = 1/512;     % Step size = 1/128
RCTFLAF_par.mu_a = 1/512;       % Step size = 1/64
P_t = 3;
RCTFLAF_par.Q = 2*P_t+1;

% Wiener TFLAF
WTFLAF_par.mu_w = 0.00390625;       % Step size
WTFLAF_par.mu_a = 0.00390625;       % Step size
P_t = 3;
WTFLAF_par.Q = 2*P_t+1;

% SFLAF parameters
SFLAF_par.mu = 0.00035;     % Step size = 1/128
P_t = 4;
SFLAF_par.L = L;
SFLAF_par.Q_t = 2*P_t+1;
SFLAF_par.expL = L*SFLAF_par.Q_t;
SFLAF_par.mu_lin = 0.001;

% PFLAF parameters
PFLAF_par.mu = 1;     % Step size = 1/128
P_t = 4;
PFLAF_par.alpha = 0.1;
PFLAF_par.Q_t = 2*P_t+1;
PFLAF_par.expL = L*PFLAF_par.Q_t;

% PHBO-FLAF parameters
PHBOFLAF_par.mu_w = 0.05;     % Step size = 1/128
PHBOFLAF_par.mu_a = 0.0005;       % Step size = 1/64
P_t = 4;
PHBOFLAF_par.Q = 2*P_t+1;
PHBOFLAF_par.alpha = 0.1;

% PSFLAF parameters
PSFLAF_par.mu = 1.5;     % Step size = 1/128
P_t = 4;
PSFLAF_par.alpha = 0.1;
PSFLAF_par.L = L;
PSFLAF_par.Q_t = 2*P_t+1;
PSFLAF_par.expL = L*PSFLAF_par.Q_t;
PSFLAF_par.mu_lin = 0.001;

% PHSAFLMS parameters
PHSAFLMS_par.mu_w_HSAFLMS = 0.1;       % Step size
PHSAFLMS_par.mu_q_HSAFLMS = 0.005;       % Step size
PHSAFLMS_par.P = WSAFLMS_par.P;      
PHSAFLMS_par.Q = WSAFLMS_par.Q;
PHSAFLMS_par.ip_range = WSAFLMS_par.ip_range;
PHSAFLMS_par.delX = WSAFLMS_par.delX ;
PHSAFLMS_par.C = WSAFLMS_par.C;

% Number of samples used to estimate steady-state MSE
mse_num_samples = 1000;

% Error variables
error_mc_LMS = zeros(1,N);
error_mc_RFFKLMS = zeros(1,N);
error_mc_WSAFLMS = zeros(1,N);
error_mc_HSAFLMS = zeros(1,N);
error_mc_TFLAF = zeros(1,N);
error_mc_TFLAF_orig = zeros(1,N);
error_mc_RCTFLAF = zeros(1,N);
error_mc_WTFLAF = zeros(1,N);
error_mc_FxdLMS = zeros(1,N);
error_mc_FxdTFLAF = zeros(1,N);
error_mc_FxdRCTFLAF = zeros(1,N);
error_mc_FxdHSAFLMS = zeros(1,N);
msd_RCTFLAF = zeros(1,N);

error_mc_SFLAF = zeros(1,N);
error_mc_PFLAF = zeros(1,N);
error_mc_PHBOFLAF = zeros(1,N);
error_mc_PSFLAF = zeros(1,N);
error_mc_PHSAFLMS = zeros(1,N);

% To write I/O to files for RTL simulation
% fid_error_rtl = fopen('/home/user/pavan/projects/adaptfilt/ASIC/LMS/simulation/inputs/error_values.txt','r');
fpath = '/home/user/pavan/projects/adaptfilt/ASIC/Retimed/LogRCTFLAF_LUT/simulation/inputs/L128_RIR_v2/';
file_write = 0;

% sin16 = rtl_bin2dec_err('/home/user/pavan/projects/adaptfilt/ASIC/TFLAF/simulation/outputs/cos_theta_out.txt',15,16);
% cos16 = rtl_bin2dec_err('/home/user/pavan/projects/adaptfilt/ASIC/TFLAF/simulation/outputs/sin_theta_out.txt',15,16);

% Load Taylor sin and cos values for fixed point simulation
% sin16 = load('/home/user/pavan/phd/work/matlab/data/error_plot_nonl/sin_val.mat');
% cos16 = load('/home/user/pavan/phd/work/matlab/data/error_plot_nonl/cos_val.mat');
% 
log_lut12 = load('/home/user/pavan/phd/work/matlab/data/LUT/log_mitchell_lut_Q4_12_no_rnd.mat');
log_lut12 = [0, log_lut12.log_mitchell_lut];
% 
% antilog_lut11 = load("/home/user/pavan/phd/work/matlab/data/LUT/antilog_mitchell_lut_Q17_11.mat");
% antilog_lut11 = antilog_lut11.antilog_lut11;

antilog_lut12 = load("/home/user/pavan/phd/work/matlab/data/LUT/antilog_mitchell_lut_Q6_12.mat");
antilog_lut12 = antilog_lut12.antilog_lut12;

% LUT logsin and logcos fxd
% lut_fraclen = 7;
% sin_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/log2_sin_lut_%dquant_Q4_12.mat', lut_fraclen);
% cos_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/log2_cos_lut_%dquant_Q4_12.mat', lut_fraclen);
% sin16 = load(sin_fpath);
% cos16 = load(cos_fpath);
% sin16 = sin16.log2sin_lut7_12;
% cos16 = cos16.log2cos_lut7_12;

% lut_fraclen = 7;
% sin_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/log2_sin_lut_%dquant_Q4_12.mat', lut_fraclen);
% cos_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/log2_cos_lut_%dquant_Q4_12.mat', lut_fraclen);
% logsin16 = load(sin_fpath);
% logcos16 = load(cos_fpath);
% logsin16 = logsin16.log2sin_lut7_12;
% logsin16(1) = 0;
% logcos16 = logcos16.log2cos_lut7_12;
% logcos16(end) = 0;

log_fraclen = 7;
% log_acc_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/log2_acc_lut_%dquant_Q5_12.mat', log_fraclen);
% log_acc_lut12 = load(log_acc_fpath);
% log_acc_lut12 = log_acc_lut12.log2_acc_lut;
log_acc_lut12 = log_lut12;

% log_acc_lut12 = load('/home/user/pavan/phd/work/matlab/data/LUT/log_mitchell_lut_Q2_14.mat');
% log_acc_lut12 = [0, log_acc_lut12.log_mitchell_lut];

% LUT Phi
lut_fraclen = 7;
sin_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/sin_lut_%dquant_Q1_15.mat', lut_fraclen);
cos_fpath = sprintf('/home/user/pavan/phd/work/matlab/data/LUT/cos_lut_%dquant_Q1_15.mat', lut_fraclen);
sin16 = load(sin_fpath);
sin16 = sin16.sin_lut7_15;
cos16 = load(cos_fpath);
cos16 = cos16.cos_lut7_15;


%     global maxval;
%     global maxvar;
%     maxval = 0;

for trial=1:monte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Generate binary data
%     x = 0.25*randn(1,N);

    % File ops
    if (file_write)
        fid_d = fopen(sprintf([fpath 'd_values%i.txt'],trial),'w');
        fid_x = fopen(sprintf([fpath 'x_values%i.txt'],trial),'w');
        fid_error = fopen(sprintf([fpath 'error_matlab%i.txt'],trial),'w');
    end

    x_white = 0.25*randn(1,N);  % Nonl filter input
%     x_white = 2*rand(1,N)-1;  % Nonl filter input
%     x_white = 1*randn(1,N);      % Lin only filter input
    
    % Thresholding to limit sample range
    th = 1.3;
    for i=1:N
        if( abs(x_white(i)) > th )
%             x_white(i) = sign(x_white(i))*th;
            x_white(i) = 0.25*x_white(i);
        end
    end
    x = x_white;
    ns = nv*randn(1,N);

    RFF_par.muRFFKLMS = 0.8;
    RFF_par.D=16;              
    RFF_par.sigma = 0.8;
    RFF_par.R = mvnrnd(zeros(L,1),(1/(RFF_par.sigma^2))*eye(L),RFF_par.D);
    RFF_par.b = 2*pi*rand(RFF_par.D,1);


%     % Coloured signal generation
%     x = zeros(1,N);
%     x(1) = x_white(1);
%     theta = 0.8;
% 
%     for i=2:N
%         x(i) = sqrt(1-theta^2)*x_white(i) + theta*x(i-1);
%     end

%     d(1) = x(1) + nv*randn(1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Linear channel for Wiener type systems

% %     x_lin = conv(x,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2,zeros(1,40),0.1,zeros(1,20)]); % L = 128
%     x_lin = conv(x,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2]);   % L
% %     = 64
% %     x_lin = conv(x,[zeros(1,5),1,0.5]);
%     x_lin = x_lin(1:length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Nonlinearity for NAEC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Asymmetrical loudspeaker nonlinearity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gain = 4;
%     d_nl = zeros(1,N);

% %       Static nonlinearity
%     z_l = 1.5*(x) - 0.3*(x).^2;

% %       Dynamic nonlinearity
    z_l = 3/2*x - 3/10*x.^2 + 9/5*x.*[0,x(1:end-1)] + 1/2*x.*[0,0,x(1:end-2)] - ... 
        2/5*x.*[0,0,0,x(1:end-3)] - 3/2*[0,x(1:end-1)].*[0,0,x(1:end-2)] + ...
        9/10*x.*[0,x(1:end-1)].*[0,0,0,x(1:end-3)] + 1/2*[0,x(1:end-1)].*[0,0,x(1:end-2)].*[0,0,0,x(1:end-3)] - ... 
        1/10*[0,x(1:end-1).^2] +1/5*[0,0,x(1:end-2).^2] - 1/10.*[0,0,0,x(1:end-3).^2] + ...
        3/10*[0,0,0,x(1:end-3)].*[0,0,0,0,x(1:end-4)] - 6/5*[0,0,0,0,x(1:end-4).^2] + ...
        1/5*[0,x(1:end-1)].*[0,0,0,0,0,x(1:end-5)] + 3/10*[0,0,0,x(1:end-3)].*[0,0,0,0,0,x(1:end-5)] + ... 
        6/5*[0,0,0,0,0,x(1:end-5).^2];
    
%     z_l = 1.5*(x+0.5) - 0.3*(x+0.5).^2;   % Shifted Q point
    rho = 0.5*(z_l<=0) + 4*(z_l>0);
    d_nl = gain*(1./(1+exp(-rho.*z_l))-0.5);
    d_nl = (d_nl-1);  % Shifting back to AC level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Symmetrical soft clipping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     zeta = 0.1;
%     d_nl = zeros(1,N);
% 
%     for ind_nl = 1:length(x)
%         if ( abs(x(ind_nl)) <= zeta )
%             d_nl(ind_nl) = 2/(3*zeta) * x(ind_nl);
%         elseif ( abs(x(ind_nl)) > zeta && abs(x(ind_nl)) <= 2*zeta )
%             d_nl(ind_nl) = sign(x(ind_nl)) * (3-(2-abs(x(ind_nl))/zeta)^2)/3;
%         else
%             d_nl(ind_nl) = sign(x(ind_nl));
%         end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loudpsk + soft clip tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %     Asymmetrical loudspeaker nonlinearity
% % 
%     gain = 2;
%     d_nl = zeros(1,N);
% 
%     z_l = 1.5*(x(1:N/2-1)) - 0.3*(x(1:N/2-1)).^2;
% %     z_l = 3/2*x - 3/10*x.^2 + 9/5*x.*[0,x(1:end-1)] + 1/2*x.*[0,0,x(1:end-2)] - ... 
% %         2/5*x.*[0,0,0,x(1:end-3)] - 3/2*[0,x(1:end-1)].*[0,0,x(1:end-2)] + ...
% %         9/10*x.*[0,x(1:end-1)].*[0,0,0,x(1:end-3)] + 1/2*[0,x(1:end-1)].*[0,0,x(1:end-2)].*[0,0,0,x(1:end-3)] - ... 
% %         1/10*[0,x(1:end-1).^2] +1/5*[0,0,x(1:end-2).^2] - 1/10.*[0,0,0,x(1:end-3).^2] + ...
% %         3/10*[0,0,0,x(1:end-3)].*[0,0,0,0,x(1:end-4)] - 6/5*[0,0,0,0,x(1:end-4).^2] + ...
% %         1/5*[0,x(1:end-1)].*[0,0,0,0,0,x(1:end-5)] + 3/10*[0,0,0,x(1:end-3)].*[0,0,0,0,0,x(1:end-5)] + ... 
% %         6/5*[0,0,0,0,0,x(1:end-5).^2];
%     
% %     z_l = 1.5*(x+0.5) - 0.3*(x+0.5).^2;   % Shifted Q point
%     rho = 0.5*(z_l<=0) + 4*(z_l>0);
%     d_nl(1:N/2-1) = gain*(1./(1+exp(-rho.*z_l))-0.5);
% %     d_nl = (d_nl-1);  % Shifting back to AC level
% 
% %     Symmetrical soft clipping
% 
%     zeta = 0.25;
% 
%     for ind_nl = N/2:length(x)
%         if ( abs(x(ind_nl)) <= zeta )
%             d_nl(ind_nl) = 2/(3*zeta) * x(ind_nl);
%         elseif ( abs(x(ind_nl)) > zeta && abs(x(ind_nl)) <= 2*zeta )
%             d_nl(ind_nl) = sign(x(ind_nl)) * (3-(2-abs(x(ind_nl))/zeta)^2)/3;
%         else
%             d_nl(ind_nl) = sign(x(ind_nl));
%         end
%     end

%     % Linear part (Hamm)
%     d = conv(d_nl,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2,zeros(1,20),0.07]);
% %     d = conv(d_nl,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2]);
%     d = d(1:length(d_nl)) + ns;

    % Wiener
%     d = d_nl + ns;

    % LMS sys
%     d = conv(x,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2,zeros(1,15),0.1]);
%     d = d(1:length(x)) + ns;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upper half buffer loudspk. Lower half buffer soft clip (Type 1 memory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     rir = fast_ISM_RoomResp(8000,[0.75 0.5 0.6 0.7 0.65 0.6],'t60', 0.06, [1 1 1.5], [2 2 1.5], [4 4 3]);
% % %     d = conv(d_nl,rir(1:512)./0.05599);         % 512 tap system
% %     d = conv(d_nl,rir(31:31+100)./0.05599);         % 128 tap system
% 
%     zeta = 0.25;
%     d_nl2 = zeros(1,length(x));
%     
%     for ind_nl = 1:length(x)
%         if ( abs(x(ind_nl)) <= zeta )
%             d_nl2(ind_nl) = 2/(3*zeta) * x(ind_nl);
%         elseif ( abs(x(ind_nl)) > zeta && abs(x(ind_nl)) <= 2*zeta )
%             d_nl2(ind_nl) = sign(x(ind_nl)) * (3-(2-abs(x(ind_nl))/zeta)^2)/3;
%         else
%             d_nl2(ind_nl) = sign(x(ind_nl));
%         end
%     end
% 
%     x_vec_gen = zeros(1,101);
%     h = rir(31:31+100)./0.05599;
%     for i = 1:N
% 
%         if(i > 51)
%             x_vec_gen = [d_nl(i),x_vec_gen(1:50), d_nl2(i-51),x_vec_gen(52:end-1)];
%         else
%             x_vec_gen = [d_nl(i),x_vec_gen(1:end-1)];
%         end
% 
%         d(i) = x_vec_gen*h;
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ceil system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d_nl = 1/8 + ceil(2*x-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple abs, square, cubic systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %     Modulus y = |x|
%     d_nl = abs(x);

%     % Parabola y = x^2
%     d_nl = x.^2;
% 
%     % Cubic y = x^3
%     d_nl = x.^3;
% 
%       Discontinuous
%     for i_s=1:N
%         
%         if abs(x(i_s)) < 0.8
%             d_nl(i_s) = x(i_s).^2;
%         else
%             d_nl(i_s) = -abs(x(i_s)).^0.5;
%         end
% 
%     end

%   cuberoot(x)

%     d_nl = nthroot(x,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AEFLN cos memoryless system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d_nl = exp(0.5*x).*(sin(pi*x) + 0.3*sin(3*pi*x) + 0.1*sin(5*pi*x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Room impulse response using RIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rir = fast_ISM_RoomResp(8000,[0.75 0.5 0.6 0.7 0.65 0.6],'t60', 0.06, [1 1 1.5], [2 2 1.5], [4 4 3]);
% %     rir(abs(rir) < 0.07*0.05599) = 0;
% %     d = conv(d_nl,rir(1:512)./0.05599);         % 512 tap system
%     d = conv(d_nl,rir(31:31+100)./0.05599);         % 128 tap system
    
%     H = load('/home/user/pavan/phd/work/colleagues/subrahmanyam_sir/Impluse_Response.mat');
%     h = H.Impulse_Resp;
%     d = conv(d_nl, h);

%     d = conv(d_nl,[zeros(1,5),1,zeros(1,15),0.5,zeros(1,10),0.2,zeros(1,20),0.07]);

    d = conv(d_nl,[1,0.5,0.2,0.07]);
%     d = conv(d_nl,[1,0,0.5,0,0.2,0,0.07]);

    d = d(1:length(d_nl)) + ns;
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rir = fast_ISM_RoomResp(8000,[0.75 0.5 0.6 0.7 0.65 0.6],'t60', 0.06, [1 1 1.5], [2 2 1.5], [4 4 3]);
%     rir2 = [zeros(20,1);rir];
% 
%     d_tmp = conv(d_nl,rir(1:512)./0.05599);         % 512 tap system
% %         d = conv(d_nl,rir(31:31+100)./0.05599);         % 64 tap system
%     d1 = d_tmp(1:length(d_nl)/2) + ns(1:length(d_nl)/2);
% 
%     d_tmp = conv(d_nl,rir2(1:512)./0.05599);         % 512 tap system
%     d_tmp = conv(d_nl,rir2(31:31+100)./0.05599);         % 64 tap system
%     d2 = d_tmp(1:length(d_nl)/2) + ns(length(d_nl)/2+1:end);
%     d = [d1, d2];

%     d = zeros(1,length(x));

    d_fi = d .* 2^(fraclen);
    d_fi = round(d_fi);
    if(overflow_det(d,fraclen,wrdlen))
        d
    end

    x_fi = x .* 2^(fraclen);
    x_fi = round(x_fi);
    if(overflow_det(x,fraclen,wrdlen))
        x
    end

    if (file_write)
        for f_ind=1:N
            fprintf(fid_d, "%s ", dec2hex(d_fi(f_ind),4));
            fprintf(fid_x, "%s ", dec2hex(x_fi(f_ind),4));
        end
    end

    % Output variables
    y_LMS = zeros(1,N);
    y_RFFKLMS = zeros(1,N);
    y_WSAFLMS = zeros(1,N);
    y_HSAFLMS = zeros(1,N);
    y_TFLAF = zeros(1,N);
    y_SFLAF = zeros(1,N);
    y_PFLAF = zeros(1,N);
    y_PHBOFLAF = zeros(1,N);
    y_PSFLAF = zeros(1,N);
    y_PHSAFLMS = zeros(1,N);

    % Error variables
    e_LMS = zeros(1,N);
    e_RFFKLMS = zeros(1,N);
    e_WSAFLMS = zeros(1,N);
    e_HSAFLMS = zeros(1,N);
    e_FxdHSAFLMS = zeros(1,N);
    e_TFLAF = zeros(1,N);
    e_TFLAF_orig = zeros(1,N);
    e_RCTFLAF = zeros(1,N);
    e_WTFLAF = zeros(1,N);
    
    e_SFLAF = zeros(1,N);
    e_PFLAF = zeros(1,N);
    e_PHBOFLAF = zeros(1,N);
    e_PSFLAF = zeros(1,N);
    e_PHSAFLMS = zeros(1,N);

    e_FxdLMS = zeros(1,N);
    e_FxdTFLAF = zeros(1,N);
    e_FxdRCTFLAF = zeros(1,N);

    e_LMS(1) = d(1);
    e_RFFKLMS(1) = d(1);
    e_WSAFLMS(1) = d(1);
    e_HSAFLMS(1) = d(1);
    e_FxdHSAFLMS(1) = d(1);
    e_TFLAF(1) = d(1);
    e_TFLAF_orig(1) = d(1);
    e_RCTFLAF(1) = d(1);
    e_WTFLAF(1) = d(1);

    e_SFLAF(1) = d(1);
    e_PFLAF(1) = d(1);
    e_PHBOFLAF(1) = d(1);
    e_PSFLAF(1) = d(1);
    e_PHSAFLMS(1) = d(1);

    e_FxdLMS(1) = nearest(d(1)*(2^(fraclen))).*(2^(-fraclen));
    e_FxdTFLAF(1) = nearest(d(1)*(2^(fraclen))).*(2^(-fraclen));
    e_FxdRCTFLAF(1) = nearest(d(1)*(2^(fraclen))).*(2^(-fraclen));

    msd_RCTFLAF = zeros(1,N);

    % Weight initialize
    w_LMS = zeros(1,L+1);
    w_FxdLMS = zeros(1,L);
    
    w_RFFKLMS = zeros(RFF_par.D+1, 1);
    
    w_WSAFLMS = zeros(1,L);
    w_WSAFLMS(1) = 1;
    q_wsaf = (-WSAFLMS_par.ip_range/2:WSAFLMS_par.delX:WSAFLMS_par.ip_range/2+5*WSAFLMS_par.delX)';
    
    w_HSAFLMS = zeros(1,L+1);
    q_hsaf = (-HSAFLMS_par.ip_range/2-HSAFLMS_par.delX:HSAFLMS_par.delX:HSAFLMS_par.ip_range/2+4*HSAFLMS_par.delX)';
    w_FxdHSAFLMS = w_HSAFLMS;
    q_hsafFxd = q_hsaf;

    w_TFLAF = zeros(1,TFLAF_par.expL+1);
    w_FxdTFLAF = zeros(1,TFLAF_par.expL+1);
    w_TFLAF_orig = zeros(1,TFLAF_par.expL+1);
    g_vec = zeros(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    g_vec_fxd = zeros(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    g_vec_orig = zeros(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    log_g_vec = zeros(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    log_g_vec_sign = zeros(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    log_g_vec_zero = ones(1,TFLAF_par.expL+del*TFLAF_par.Q_t);
    
    w_RCTFLAF = zeros(1,L+1);
    a_RCTFLAF = zeros(1,RCTFLAF_par.Q);
    a_RCTFLAF(1) = 1;
    s_vec_RCTFLAF = zeros(1,L+del);
    g_vec_RCTFLAF = zeros(L+del,RCTFLAF_par.Q);

    w_FxdRCTFLAF = zeros(1,L+1);
    a_FxdRCTFLAF = zeros(1,RCTFLAF_par.Q);
    a_FxdRCTFLAF(1) = 1;
    del_a_upd = 2;
    a_FxdRCTFLAF_del = zeros(del_a_upd+1,RCTFLAF_par.Q);
    for var_del = 1:del_a_upd
        a_FxdRCTFLAF_del(var_del,:) = a_FxdRCTFLAF;
    end
    log_s_vec_fxd_RCTFLAF = zeros(1,L+del);
    log_g_vec_fxd_RCTFLAF = zeros(L+del,RCTFLAF_par.Q);
    log_s_vec_sign_fxd = zeros(1,L+del);
    log_s_vec_zero_fxd = ones(1,L+del);
    log_g_vec_sign_fxd = zeros(L+del,RCTFLAF_par.Q);
    log_g_vec_zero_fxd = ones(L+del,RCTFLAF_par.Q);

    w_WTFLAF = zeros(1,L+1);
    w_WTFLAF(2) = 1;
    a_WTFLAF = zeros(1,WTFLAF_par.Q);
    a_WTFLAF(1) = 1;
    s_vec_WTFLAF = zeros(1,L);

    % HSAF init
    s_vec = zeros(1,L+del);
    U_vec = zeros(HSAFLMS_par.P+1,L+del);
    U_vec_C_fxd = zeros(HSAFLMS_par.P+1,L+del);
    j_vec = zeros(1,L+del);
    s_vec_fxd = zeros(1,L+del);
    U_vec_fxd = zeros(HSAFLMS_par.P+1,L+del);
    j_vec_fxd = zeros(1,L+del);

    % PHSAF init
    s_vec_PHSAF = zeros(1,L+del);
    U_vec_PHSAF = zeros(PHSAFLMS_par.P+1,L+del);
    j_vec_PHSAF = zeros(1,L+del);

%     [e_HSAFLMS(1), w_HSAFLMS, q_hsaf, s_vec, U_vec, j_vec] = HSAFLMS_filter(x(1), s_vec, U_vec, j_vec, d(1), w_HSAFLMS, q_hsaf, HSAFLMS_par);
    

    w_SFLAF = zeros(1,SFLAF_par.expL+1);
    w_SFLAF_lin = zeros(1,SFLAF_par.L+1);
    g_vec_SFLAF = zeros(1,SFLAF_par.expL+del*SFLAF_par.Q_t);
    w_PFLAF = zeros(1,PFLAF_par.expL+1);
    g_vec_PFLAF = zeros(1,PFLAF_par.expL+del*PFLAF_par.Q_t);
    w_PHBOFLAF = zeros(1,L+1);
    a_PHBOFLAF = zeros(1,PHBOFLAF_par.Q);
    a_PHBOFLAF(1) = 1;
    s_vec_PHBOFLAF = zeros(1,L+del);
    g_vec_PHBOFLAF = zeros(L+del,PHBOFLAF_par.Q);
    w_PSFLAF = zeros(1,PSFLAF_par.expL+1);
    w_PSFLAF_lin = zeros(1,PSFLAF_par.L+1);
    g_vec_PSFLAF = zeros(1,PSFLAF_par.expL+del*PSFLAF_par.Q_t);

    w_PHSAFLMS = zeros(1,L+1);
    q_hsaf_PHSAF = (-PHSAFLMS_par.ip_range/2-PHSAFLMS_par.delX:PHSAFLMS_par.delX:PHSAFLMS_par.ip_range/2+4*PHSAFLMS_par.delX)';

    % Channel related
    %x_vec = zeros(1,L);
    x_vec = zeros(1,L+del);     % Delayed version
    %x_vec(1) = x(1);

    z_op = zeros(1,N);
    z_op(1) = x(1);
    j_track = zeros(1,N);
    delw_track = zeros(N,TFLAF_par.expL+1);  
    q_n_sum_track = zeros(1,N);
    s_vec2_RC_track = zeros(1,N);
    Gw_RC_track = zeros(1,N);

    s_vec2_HSAF_track = zeros(1,N);
    ucw_HSAF_track = zeros(1,N);

    % Start MC experiment
    for i=1:N

        % Non-linear channels
        % RFF paper example
%         d(i) = (0.9  - ( 0.5* exp(-(x(i)^2)))*x(i-1) + 0.2*sqrt((sin(pi*x((i-2)))^2))) + nv*randn(1,1); 

        % Example 1
%         z_op(i) = 1.1 * exp(-abs(z_op(i-1))) + x(i); 
%         d(i) = z_op(i).^2 + nv*randn(1,1); 
        % Example 2 
%         if(i > 2 && i < 10000)
%             d(i) = 0.1*sin(d(i-1)*pi) + 0.4*cos(d(i-2)*pi) + nv*randn(1,1);
%         elseif( i >= 10000)
%             d(i) = 0.9*sin(d(i-1)*pi) + 0.15*cos(d(i-2)*pi) + nv*randn(1,1);
%         end
          % Example 3
%         d(i) = d(i-1)/(1+d(i-1)^2)+x(i-1)^3 + nv*randn(1,1);

        % AEFLN system (Type 1 memory system)
%         d(i) = 0.6*sin(pi*x(i))^3 - 2/(x(i).^3+2) - 0.1*cos(4*pi*x(i-4)) + 1.125 + nv*randn(1,1);
%         d(i) = 0.6*sin(pi*x(i))^3 + 0.2*cos(2*pi*x(i-2))^2 - 0.1*cos(4*pi*x(i-4)) + 1.125 + nv*randn(1,1);

        % Input tap
        x_vec = [x(i), x_vec(1:end-1)];

%         % Apply to filters
% % 
%         [e_LMS(i)] = LMS_filter(x_vec(1:L), d(i), w_LMS);
% %         [e_FxdLMS(i)] = LMS_filter_fxd(x_vec(1:L), d(i), w_FxdLMS, fraclen, wrdlen);
%         [e_TFLAF(i), g_vec] = TFLAF_filter(x(i), d(i), w_TFLAF, g_vec, TFLAF_par);
%         [e_FxdTFLAF(i), log_g_vec, log_g_vec_sign, log_g_vec_zero] = LogTFLAF_filter_fxd(x(i), d(i), w_FxdTFLAF, log_g_vec, log_g_vec_zero, log_g_vec_sign, TFLAF_par, fraclen, lut_fraclen, wrdlen, logsin16, logcos16, log_lut12, antilog_lut12, log_acc_lut12);
%         [e_FxdTFLAF(i), log_g_vec, log_g_vec_sign, log_g_vec_zero] = LogTFLAF_filter_fxd(x(i), d(i), w_FxdTFLAF, log_g_vec, log_g_vec_zero, log_g_vec_sign, TFLAF_par, fraclen, lut_fraclen, wrdlen, sin16, cos16, log_lut12, antilog_lut12, log_acc_lut12);
%           [e_FxdTFLAF(i), g_vec_fxd] = TFLAF_filter_fxd(x(i), d(i), w_FxdTFLAF, g_vec_fxd, TFLAF_par, fraclen, wrdlen, sin16, cos16);
%         [e_RCTFLAF(i), s_vec_RCTFLAF, g_vec_RCTFLAF] = RCTFLAF_filter(x(i), d(i), w_RCTFLAF, a_RCTFLAF, s_vec_RCTFLAF, g_vec_RCTFLAF, RCTFLAF_par);
%         s_vec2_RC_track(i) = s_vec_RCTFLAF*s_vec_RCTFLAF';
%         Gw_RC_track(i)     = norm(w_RCTFLAF(2:end)*g_vec_RCTFLAF)^2;      
%         [e_FxdRCTFLAF(i), log_s_vec_fxd_RCTFLAF, log_s_vec_sign_fxd, log_s_vec_zero_fxd, log_g_vec_fxd_RCTFLAF, log_g_vec_sign_fxd, log_g_vec_zero_fxd] = LogRCTFLAF_filter_fxd(x(i), d(i), w_FxdRCTFLAF, a_FxdRCTFLAF, log_s_vec_fxd_RCTFLAF, log_s_vec_sign_fxd, log_s_vec_zero_fxd, log_g_vec_fxd_RCTFLAF, log_g_vec_sign_fxd, log_g_vec_zero_fxd, RCTFLAF_par, fraclen, wrdlen, logsin16, logcos16, antilog_lut12, log_lut12, lut_fraclen, log_acc_lut12, fraclen);
%         [e_FxdRCTFLAF(i), log_s_vec_fxd_RCTFLAF, log_s_vec_sign_fxd, log_s_vec_zero_fxd, log_g_vec_fxd_RCTFLAF, log_g_vec_sign_fxd, log_g_vec_zero_fxd] = LogRCTFLAF_filter_fxd(x(i), d(i), w_FxdRCTFLAF, a_FxdRCTFLAF, log_s_vec_fxd_RCTFLAF, log_s_vec_sign_fxd, log_s_vec_zero_fxd, log_g_vec_fxd_RCTFLAF, log_g_vec_sign_fxd, log_g_vec_zero_fxd, RCTFLAF_par, fraclen, wrdlen, sin16, cos16, antilog_lut12, log_lut12, lut_fraclen, log_acc_lut12, log_fraclen);
        [e_HSAFLMS(i), s_vec, U_vec, j_vec] = HSAFLMS_filter(x(i), s_vec, U_vec, j_vec, d(i), w_HSAFLMS, q_hsaf, HSAFLMS_par);
%         s_vec2_HSAF_track(i) = s_vec*s_vec';
%         ucw_HSAF_track(i)    = norm(HSAFLMS_par.C'*U_vec*w_HSAFLMS(2:end)')^2;

        [e_FxdHSAFLMS(i), s_vec_fxd, U_vec_fxd, U_vec_C_fxd, j_vec_fxd] = HSAFLMS_filter_fxd(x(i), s_vec_fxd, U_vec_fxd, U_vec_C_fxd, j_vec_fxd, d(i), w_FxdHSAFLMS, q_hsafFxd, HSAFLMS_par, fraclen, wrdlen);
% % 
% %         [e_SFLAF(i), g_vec_SFLAF] = SFLAF_filter(x(i), x_vec(1:L), d(i), w_SFLAF, w_SFLAF_lin, g_vec_SFLAF, SFLAF_par);
%         [e_PFLAF(i), g_vec_PFLAF] = PFLAF_filter(x(i), d(i), w_PFLAF, g_vec_PFLAF, PFLAF_par);
%         [e_PHBOFLAF(i), s_vec_PHBOFLAF, g_vec_PHBOFLAF] = PHBOFLAF_filter(x(i), d(i), w_PHBOFLAF, a_PHBOFLAF, s_vec_PHBOFLAF, g_vec_PHBOFLAF, PHBOFLAF_par);
% %         [e_PSFLAF(i), g_vec_PSFLAF] = PSFLAF_filter(x(i), x_vec(1:L), d(i), w_PSFLAF, w_PSFLAF_lin, g_vec_PSFLAF, PSFLAF_par);
%         [e_PHSAFLMS(i), s_vec_PHSAF, U_vec_PHSAF, j_vec_PHSAF] = PHSAFLMS_filter(x(i), s_vec_PHSAF, U_vec_PHSAF, j_vec_PHSAF, d(i), w_PHSAFLMS, q_hsaf_PHSAF, PHSAFLMS_par);

        % Apply weight update

        if (i > del)
%             [w_LMS] = LMS_update(x_vec(del+1:end), w_LMS, mu_LMS, e_LMS(i-del));
%             [w_FxdLMS] = LMS_update_fxd(x_vec(del+1:end), w_FxdLMS, mu_LMS, e_FxdLMS(i-del), fraclen, wrdlen);
%             [w_TFLAF] = TFLAF_update(w_TFLAF, e_TFLAF(i-del), g_vec(del*TFLAF_par.Q_t+1:end), TFLAF_par);
%             [w_FxdTFLAF] = LogTFLAF_update_fxd(w_FxdTFLAF, e_FxdTFLAF(i-del), log_g_vec(del*TFLAF_par.Q_t+1:end), log_g_vec_zero(del*TFLAF_par.Q_t+1:end), log_g_vec_sign(del*TFLAF_par.Q_t+1:end), TFLAF_par, fraclen, wrdlen, log_lut12, antilog_lut12);

%             [w_FxdTFLAF,delw_track(i,:)] = TFLAF_update_fxd(w_FxdTFLAF,e_FxdTFLAF(i-del),g_vec_fxd(del*TFLAF_par.Q_t+1:end), TFLAF_par, fraclen, wrdlen);
% %             w_old = a_RCTFLAF;
%             [w_RCTFLAF, a_RCTFLAF] = RCTFLAF_update(w_RCTFLAF, a_RCTFLAF, s_vec_RCTFLAF(del+1:end), g_vec_RCTFLAF(del+1:end,:), e_RCTFLAF(i-del), RCTFLAF_par);
% % %             w_cmp(i,:) = w_old - a_RCTFLAF;                                       
%             [w_FxdRCTFLAF, a_FxdRCTFLAF] = LogRCTFLAF_update_fxd(w_FxdRCTFLAF, a_FxdRCTFLAF, log_s_vec_fxd_RCTFLAF(del+1:end), log_s_vec_sign_fxd(del+1:end), log_s_vec_zero_fxd(del+1:end), log_g_vec_fxd_RCTFLAF(del+1:end,:), log_g_vec_sign_fxd(del+1:end,:), log_g_vec_zero_fxd(del+1:end,:), e_FxdRCTFLAF(i-del), RCTFLAF_par, fraclen, wrdlen, log_lut12, antilog_lut12, log_acc_lut12, fraclen);
%             [w_FxdRCTFLAF, a_FxdRCTFLAF] = RCTFLAF_update_fxd_log_vmm(w_FxdRCTFLAF, a_FxdRCTFLAF, s_vec_fxd_RCTFLAF(del+1:end), g_vec_fxd_RCTFLAF(del+1:end,:), e_FxdRCTFLAF(i-del), RCTFLAF_par, fraclen, wrdlen, log_lut12, antilog_lut11);
            [w_HSAFLMS, q_hsaf] = HSAFLMS_update(s_vec(del+1:end), U_vec(:,del+1:end), j_vec(del+1:end), e_HSAFLMS(i-del), w_HSAFLMS, q_hsaf, HSAFLMS_par);
              [w_FxdHSAFLMS, q_hsafFxd] = HSAFLMS_update_fxd(s_vec_fxd(del+1:end), U_vec_fxd(:,del+1:end), U_vec_C_fxd(:,del+1:end), j_vec_fxd(del+1:end), e_FxdHSAFLMS(i-del), w_FxdHSAFLMS, q_hsafFxd, HSAFLMS_par, fraclen, wrdlen);
% %               [w_FxdHSAFLMS, q_hsafFxd] = HSAFLMS_update(s_vec_fxd(del+1:end), U_vec_fxd(:,del+1:end), j_vec_fxd(del+1:end), e_FxdHSAFLMS(i-del), w_FxdHSAFLMS, q_hsafFxd, HSAFLMS_par);
% % 
% %             [w_SFLAF, w_SFLAF_lin] = SFLAF_update(x_vec(del+1:end), w_SFLAF, w_SFLAF_lin, e_SFLAF(i-del), g_vec_SFLAF(del*SFLAF_par.Q_t+1:end), SFLAF_par);
%             [w_PFLAF] = PFLAF_update(w_PFLAF, e_PFLAF(i-del), g_vec_PFLAF(del*PFLAF_par.Q_t+1:end), PFLAF_par);
%             [w_PHBOFLAF, a_PHBOFLAF] = PHBOFLAF_update(w_PHBOFLAF, a_PHBOFLAF, s_vec_PHBOFLAF(del+1:end), g_vec_PHBOFLAF(del+1:end,:), e_PHBOFLAF(i-del), PHBOFLAF_par);
% %             [w_PSFLAF, w_PSFLAF_lin] = PSFLAF_update(x_vec(del+1:end), w_PSFLAF, w_PSFLAF_lin, e_PSFLAF(i-del), g_vec_PSFLAF(del*PSFLAF_par.Q_t+1:end), PSFLAF_par);
%             [w_PHSAFLMS, q_hsaf_PHSAF] = PHSAFLMS_update(s_vec_PHSAF(del+1:end), U_vec_PHSAF(:,del+1:end), j_vec_PHSAF(del+1:end), e_PHSAFLMS(i-del), w_PHSAFLMS, q_hsaf_PHSAF, PHSAFLMS_par);
% 
        end

%         [e_RFFKLMS(i), w_RFFKLMS] = RFFKLMS_filter(x_vec(1:L), d(i), w_RFFKLMS, RFF_par);
%         [e_WSAFLMS(i), w_WSAFLMS, q_wsaf] = WSAFLMS_filter(x_vec(1:L), d(i), w_WSAFLMS, q_wsaf, WSAFLMS_par);
%             if i > 1
%             [e_HSAFLMS(i), w_HSAFLMS, q_hsaf, s_vec, U_vec, j_vec] = HSAFLMS_filter(x(i), s_vec, U_vec, j_vec, d(i), w_HSAFLMS, q_hsaf, HSAFLMS_par);
%             end
%         [e_WTFLAF(i), w_WTFLAF, a_WTFLAF] = WTFLAF_filter(x_vec(1:L), d(i), w_WTFLAF, a_WTFLAF, WTFLAF_par);

%         j_track(i) = j_vec(1);

        % RC-TFLAF MSD calculation

%         msd_RCTFLAF(i) = norm(rir(1:512)./0.05599 - w_RCTFLAF(2:end))/norm(rir(1:512)./0.05599);
%         wb_track(i) = w_RCTFLAF(1);

    end
        
    % MC error accumulation
    error_mc_LMS = error_mc_LMS + (e_LMS.^2);
    error_mc_RFFKLMS = error_mc_RFFKLMS + (e_RFFKLMS.^2);
    error_mc_WSAFLMS = error_mc_WSAFLMS + (e_WSAFLMS.^2);
    error_mc_HSAFLMS = error_mc_HSAFLMS + (e_HSAFLMS.^2);
    error_mc_TFLAF = error_mc_TFLAF + (e_TFLAF.^2);
    error_mc_TFLAF_orig = error_mc_TFLAF_orig + (e_TFLAF_orig.^2);
    error_mc_RCTFLAF = error_mc_RCTFLAF + (e_RCTFLAF.^2);
    error_mc_WTFLAF = error_mc_WTFLAF + (e_WTFLAF.^2);
    
    error_mc_SFLAF = error_mc_SFLAF + (e_SFLAF.^2);
    error_mc_PFLAF = error_mc_PFLAF + (e_PFLAF.^2);
    error_mc_PHBOFLAF = error_mc_PHBOFLAF + (e_PHBOFLAF.^2);
    error_mc_PSFLAF = error_mc_PSFLAF + (e_PSFLAF.^2);
    error_mc_PHSAFLMS = error_mc_PHSAFLMS + (e_PHSAFLMS.^2);

    error_mc_FxdLMS = error_mc_FxdLMS + (e_FxdLMS.^2);
    error_mc_FxdTFLAF = error_mc_FxdTFLAF + (e_FxdTFLAF.^2);
    error_mc_FxdRCTFLAF = error_mc_FxdRCTFLAF + (e_FxdRCTFLAF.^2);
    error_mc_FxdHSAFLMS = error_mc_FxdHSAFLMS + (e_FxdHSAFLMS.^2);

    e_FxdLMS_fi = e_FxdLMS .* 2^(fraclen);
    e_FxdLMS_fi = nearest(e_FxdLMS_fi);

    e_FxdTFLAF_fi = e_FxdTFLAF .* 2^(fraclen);
    e_FxdTFLAF_fi = nearest(e_FxdTFLAF_fi);

    e_FxdRCTFLAF_fi = e_FxdRCTFLAF .* 2^(fraclen);
    e_FxdRCTFLAF_fi = nearest(e_FxdRCTFLAF_fi);

%     msd_mc_RCTFLAF = msd_mc_RCTFLAF + msd_RCTFLAF;

    if( file_write )
        for f_ind=1:N
            fprintf(fid_error, "%s ", dec2hex(e_FxdRCTFLAF_fi(f_ind),4));
        end

        save(sprintf([fpath 'error%i.mat'],trial), 'e_FxdRCTFLAF');

        fclose(fid_d);
        fclose(fid_x);
        fclose(fid_error);
    end

end

% Compute error over all MC experiments + 20 order averaging
% error_mc_LMS = movmean(error_mc_LMS / monte,20);
% error_mc_RFFKLMS = movmean(error_mc_RFFKLMS / monte,20);
% error_mc_WSAFLMS = movmean(error_mc_WSAFLMS / monte,20);
% error_mc_HSAFLMS = movmean(error_mc_HSAFLMS / monte,20);
% error_mc_TFLAF = movmean(error_mc_TFLAF / monte,20);
% error_mc_TFLAF_orig = movmean(error_mc_TFLAF_orig / monte,20);
% error_mc_RCTFLAF = movmean(error_mc_RCTFLAF / monte,20);
% error_mc_WTFLAF = movmean(error_mc_WTFLAF / monte,20);
% 
% error_mc_FxdLMS = movmean(error_mc_FxdLMS / monte,20);
% error_mc_FxdTFLAF = movmean(error_mc_FxdTFLAF / monte,20);
% error_mc_FxdRCTFLAF = movmean(error_mc_FxdRCTFLAF / monte,20);
% error_mc_FxdHSAFLMS = movmean(error_mc_FxdHSAFLMS / monte,20);
% 
% error_mc_SFLAF = movmean(error_mc_SFLAF / monte,20);
% error_mc_PFLAF = movmean(error_mc_PFLAF / monte,20);
% error_mc_PHBOFLAF = movmean(error_mc_PHBOFLAF / monte,20);
% error_mc_PSFLAF = movmean(error_mc_PSFLAF / monte,20);
% error_mc_PHSAFLMS = movmean(error_mc_PHSAFLMS / monte,20);

error_mc_LMS = movmean(error_mc_LMS / monte,500);
error_mc_TFLAF = movmean(error_mc_TFLAF / monte,500);
error_mc_RCTFLAF = movmean(error_mc_RCTFLAF / monte,500);
error_mc_FxdTFLAF = movmean(error_mc_FxdTFLAF / monte,500);
error_mc_FxdRCTFLAF = movmean(error_mc_FxdRCTFLAF / monte,500);
error_mc_HSAFLMS = movmean(error_mc_HSAFLMS / monte,500);
error_mc_FxdHSAFLMS = movmean(error_mc_FxdHSAFLMS / monte,500);
error_mc_RFFKLMS = movmean(error_mc_RFFKLMS / monte,500);

% msd_mc_RCTFLAF = msd_mc_RCTFLAF / monte;

fpath1 = '/home/user/pavan/';

toc

% Plot figures
figure;
% plot(10*log10(error_mc_LMS));
hold on;
% % plot(10*log10(error_mc_FxdLMS));
% plot(10*log10(error_mc_RFFKLMS));
% % % plot(10*log10(error_mc_WSAFLMS));
plot(10*log10(error_mc_HSAFLMS));
plot(10*log10(error_mc_FxdHSAFLMS));
% plot(10*log10(error_mc_TFLAF));
% plot(10*log10(error_mc_TFLAF_orig));  
% plot(10*log10(error_mc_FxdTFLAF));
% plot(10*log10(error_mc_RCTFLAF));
% plot(10*log10(error_mc_FxdRCTFLAF));
% % % plot(10*log10(error_mc_WTFLAF));

% plot(10*log10(error_mc_SFLAF));
% plot(10*log10(error_mc_PFLAF));
% plot(10*log10(error_mc_PHBOFLAF),'k');
% plot(10*log10(error_mc_PSFLAF),'b');
% plot(10*log10(error_mc_PHSAFLMS),'g');

% save(sprintf([fpath1 'error_tflaf%i.mat'],fraclen),"error_mc_FxdTFLAF")
% save(sprintf([fpath1 'error_rctflaf%i.mat'],fraclen),"error_mc_FxdRCTFLAF")

grid on
xlabel('Iteration');
ylabel('MSE (dB)');
title('MSE learning curves');
% legend("LMS", "RFF-KLMS", "HSAF", "TFLAF","RC-TFLAF");
legend("LMS","HSAF", "TFLAF","RC-TFLAF");

disp("Steady state MSE ")

steadystateMSE_RCTFLAF = mean(10*log10(error_mc_RCTFLAF(end-mse_num_samples:end)))
steadystateMSE_FxdRCTFLAF = mean(10*log10(error_mc_FxdRCTFLAF(end-mse_num_samples:end)))

steadystateMSE_TFLAF = mean(10*log10(error_mc_TFLAF(end-mse_num_samples:end)))
steadystateMSE_FxdTFLAF = mean(10*log10(error_mc_FxdTFLAF(end-mse_num_samples:end)))
steadystateMSE_HSAFLMS = mean(10*log10(error_mc_HSAFLMS(end-mse_num_samples:end)))
steadystateMSE_FxdHSAFLMS = mean(10*log10(error_mc_FxdHSAFLMS(end-mse_num_samples:end)))

steadystateMSE_RFFKLMS = mean(10*log10(error_mc_RFFKLMS(end-mse_num_samples:end)))

