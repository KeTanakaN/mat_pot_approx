clear all;
dir_smp = './DataSmp/'; % Directory to save the data of the optimal sampling points

%% Functions related to a weight / Initial value of alpha / Range of x
% liw = (logarithm of the inverse of a weight function)
% dlw = (derivative of the logarithm of a weight function)
% ddlw = (2nd derivative of the logarithm of a weight function)

% [SE]
% liw = @(x) log(cosh(2*x));
% dlw = @(x) -2*tanh(2*x);
% ddlw = @(x) -4./cosh(2*x).^2;
% init_alpha = @(N) (pi/2)*sqrt(N+1);
% adj_eps = 1; % "no adjust" version
% x = [-20:0.01:20];
% prefix = 'se'; % for filenames
% Prefix = 'SE'; % for titles of graphs

% [Gauss]
% liw = @(x) x.^2;
% dlw = @(x) -2*x;
% ddlw = @(x) -2;
% init_alpha = @(N) ((3*(pi^2)*(N+1))^(1/3))/2;
% adj_eps = 1; % "no adjust" version
% x = [-7.5:0.01:7.5];
% prefix = 'ga';    % for filenames
% Prefix = 'Gauss'; % for titles of graphs

% [DE]
liw = @(x) LogInvDE(2*x);
dlw = @(x) -pi*cosh(2*x).*tanh((pi/2)*sinh(2*x));
ddlw = @(x) -pi*((pi*cosh(2*x).^2)./(cosh((pi/2)*sinh(2*x)).^2) + 2*sinh(2*x).*tanh((pi/2)*sinh(2*x)));
init_alpha = @(N) lambertw(pi*(N+1))/2;
adj_eps = 1; % "no adjust" version
x = [-2.5:0.01:2.5];
prefix = 'de'; % for filenames
Prefix = 'DE'; % for titles of graphs

% [DE (only for test)]
% liw = @(x) exp(abs(x));
% dlw = @(x) -sign(x).*exp(abs(x));
% ddlw = @(x) -exp(abs(x));
% init_alpha = @(N) lambertw(pi^2*(N+1)/6);
% adj_eps = 0.5;
% x = [-4.5:0.01:4.5];
% prefix = 'detest'; % for filenames
% Prefix = 'DEtest'; % for titles of graphs

%% Generating the optimal sampling points
i = 1;
maxlog10bw = zeros(1,10);
for N = 10:10:100 % N = ((number of the optimal sampling points)-1)/2
    %% Generating the optimal sampling points
    M = 2^12; % number of sampling points to discretize the Fourier transforms
    alpha = SUB_alpha(dlw, ddlw, init_alpha(N), N, M, adj_eps);
    [a,tmp_N] = SUB_gen_opt_sample(liw, dlw, alpha, M, adj_eps);

    %% Output and plot of the sampling points
    filename = strcat(dir_smp, prefix, '_N_', num2str(N), '.txt');
    % dlmwrite(filename,a,'precision',15);
    
    plot([-tmp_N:tmp_N], a, 'o-', 'LineWidth', 2);
    set(gca,'FontName','Times','FontSize',10,'FontWeight','bold');
    xlabel('n');
    ylabel('Sampling points');
    title(['Sampling points for the ', Prefix ,' weight (N = ',num2str(N),')']);
    grid on;
    
    pause;
    
    %% Output and plot of the transformed Blaschke product with weight
    lx = length(x);
    la = length(a);
    SumLTW = zeros(1,lx);
    for n=1:la
        SumLTW = SumLTW + log(abs(tanh(x-a(n))));
    end
    LIW = zeros(1,lx);
    for l=1:lx
        LIW(l) = liw(x(l));
    end
    SumLTW = SumLTW - LIW;
    
    maxlog10bw(i) = max(SumLTW/log(10));
    i = i+1;
    
    plot(x, SumLTW/log(10),'LineWidth', 2);
    set(gca,'FontName','Times','FontSize',12,'FontWeight','bold');
    xlabel('x');
    ylabel('Function values');
    title(['Discrete weighted potential with the ', Prefix ,' weight (N = ',num2str(N),')']);
    grid on;
    
    pause;
end

%% Output of the maximum values of 'log10(transformed Blaschke product with weight)'
% filename = strcat(dir_smp, prefix, '_MaxLog10Bw.txt');
% dlmwrite(filename, maxlog10bw,'precision',15);
