clear all;
dir_smp = './DataSmp/'; % Directory to save the data of the optimal sampling points
dir_err = './DataErr/'; % Directory to save the data of the errors 

%% This program uses ADVANPIX Multiprecision Computing Toolbox
% http://www.advanpix.com/documentation/users-manual/
addpath('C:\Users\ketanaka\Documents\Multiprecision Computing Toolbox');
mp.Digits(90);

%% Functions etc.
% [SE]
% f = @(x) cosh(2*x).^(-1);
% liw = @(x) log(cosh(2*x));
% fwnorm = 1;
% x = mp([-20:0.04:20]); % MULTI
% prefix = 'se'; % for filenames

% [Gauss]
% f = @(x) (x.^2).*exp(-x.^2)./((mp('pi')^2)/16+x.^2); % MULTI
% liw = @(x) x.^2;
% fwnorm = 1;
% x = mp([-10:0.02:10]); % MULTI
% prefix = 'ga'; % for filenames

% [DE]
f = @(x) cosh((mp('pi')/2)*sinh(2*x)).^(-1); % MULTI
liw  = @(x) log(cosh((mp('pi')/2)*sinh(2*x))); % MULTI
fwnorm = 1;
x = mp([-2.5:0.005:2.5]); % MULTI
prefix = 'de'; % for filenames

%% Approximations
i = 1;
err_SE  = zeros(1,10);
err_DE  = zeros(1,10);
err_Gan = zeros(1,10);
err_opt = zeros(1,10);
for N = 10:10:100
    tic;
    %% SE-Sinc formula
%     h = mp('pi')/(2*sqrt(2*N)); % MULTI
%     h = (mp('pi')/(2*N))^(2/3); % MULTI
%     appr_SE = SUB_approx_Sinc_MP(x, N, h, f);

    %% DE-Sinc formula
%     h = log(2*mp('pi')*N)/(2*N); % MULTI
%     appr_DE = SUB_approx_Sinc_MP(x, N, h, f);
    
    %% Ganelius formula
%     samp_Gan = mp(SUB_gen_Gan_sample(N)); % MULTI
%     appr_Gan = SUB_approx_Blaschke_MP(x, samp_Gan, f, liw); % MULTI
    
    %% Optimal formula 
    filename = strcat(dir_smp, prefix, '_N_', num2str(N), '.txt');
    samp_opt = mp(dlmread(filename)); % MULTI, [SE: se_N_xx.txt, Gaussian: ga_N_xx.txt, DE: de_N_xx.txt]
    appr_opt = SUB_approx_Blaschke_MP(x, samp_opt, f, liw); % MULTI

    %% Errors
%     err_SE(i)  = max(log10(abs(appr_SE-f(x))/fwnorm));
%     err_DE(i)  = max(log10(abs(appr_DE-f(x))/fwnorm));
%     err_Gan(i) = max(log10(abs(appr_Gan-f(x))/fwnorm));
    err_opt(i) = max(log10(abs(appr_opt-f(x))/fwnorm));
    i = i+1;

    hold off;
%     plot(x, log10(abs(appr_SE-f(x))/fwnorm), '-r', 'LineWidth', 2);
%     hold on;
%     plot(x, log10(abs(appr_DE-f(x))/fwnorm), '-b', 'LineWidth', 2);
%     plot(x, log10(abs(appr_Gan-f(x))/fwnorm), '-k', 'LineWidth', 2);
%     plot(x, log10(abs(appr_opt-f(x))/fwnorm), '-g', 'LineWidth', 2);
%     plot(x, log10(abs(appr_opt-f(x))/fwnorm) - liw(x)/log(10));
%     legend('SE-Sinc', 'DE-Sinc', 'Ganelius', 'Optimal');
    grid on;
    set(gca,'FontName','Times','FontSize',15,'FontWeight','bold');
    xlabel('x');
    ylabel('log_{10}(Absolute Error)');
%     xlim([-20,20]);
%     ylim([-25,-1]);    
    
    toc;
    
%     pause;
end

%% Output of the errors
% filename = strcat(dir_err, prefix, '_err_SE.txt');
% dlmwrite(filename, err_SE);
% filename = strcat(dir_err, prefix, '_err_DE.txt');
% dlmwrite(filename, err_DE);
% filename = strcat(dir_err, prefix, '_err_Gan.txt');
% dlmwrite(filename, err_Gan);
filename = strcat(dir_err, prefix, '_err_opt.txt');
dlmwrite(filename, err_opt);
