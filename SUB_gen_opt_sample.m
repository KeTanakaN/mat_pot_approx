%%%%% [Generator of the optimal sampling points] %%%%%
% liw = (logarithm of the inverse of a weight function)
% dlw = (derivative of the logarithm of a weight function)
% alpha = (parameter determining the support [-alpha, alpha] of the optimal measure)
% M = (number of sampling points to discretize the Fourier transforms)

function [a,N] = SUB_gen_opt_sample(liw, dlw, alpha, M, adj_eps)
    %% Setting the grids for the phisical (x) and frequency (w) spaces
    h1= alpha/M;
    x = [-alpha+h1:h1:alpha];
    h2= pi/alpha;
    w = [-(pi*M)/alpha+h2:h2:(pi*M)/alpha];

    %% Fourier transform of the (approximately) extended log-weight
    lwa = liw(x) - liw(alpha);
    wsc = dlw(alpha)*SC(alpha, w);
    wsc(M) = dlw(alpha)*SC(alpha, 0);
    ftw = h1 * FFFT(lwa, h1*h2/(2*pi), M) + wsc * adj_eps;

    %% Fourier transform of the logarithm of tanh
    flt = FLT(w);
    flt(M) = FLT(0);

    %% Solution b' by inverse Fourier transform
    bp = (h2/(2*pi)) * FFFT(ftw./flt, h1*h2/(2*pi), M);
    
    %% Indefinite integral of b'
    xarrb = [-alpha:h1:alpha];
    b = zeros(1,2*M+1);
    b(M+1)=0;
    for n=1:M
        b(M+n+1) = b(M+n) + bp(M+n)*h1;
        b(M-n+1) = -b(M+n+1);
    end
    b = real(b);
    
    %% Generating sampling points
    N = floor(max(b));
    xeval = [-N:N];
    a = interp1(b,xarrb,xeval,'pchip'); % a(-N),a(-N+1),...,a(N-1),a(N) are sampling points
end
