%%%%% [Solver of the equation of alpha] %%%%%
% dlw = (derivative of the logarithm of a weight function)
% ddlw = (2nd derivative of the logarithm of a weight function)
% init_alpha = (initial value of alpha)
% N = (parameter for the optimal sampling points)
%     (2N+1 is the number of the optimal sampling points)
% M = (number of sampling points to discretize the integral)

function alpha = SUB_alpha(dlw, ddlw, init_alpha, N, M, adj_eps)
    alpha = init_alpha;

    %% Newton method
    eps = 1e-10;
    crtn = 10 * eps;
    
    while crtn > eps
        % Setting the grids for the phisical (x) and frequency (w) spaces
        h1= alpha/M;
        x = [-alpha+h1:h1:alpha];
        
        % 1 step
        % func = -(2/pi^2)*h1*(x*dlw(x)') - ((2*alpha+1)/pi^2)*dlw(alpha) * adj_eps - (N+1);
        func = -(2/pi^2)*h1*(x*dlw(x)') - ((2*alpha+1)/pi^2)*dlw(alpha) * adj_eps - (N+adj_eps);
        % diff = - ((2*alpha+1)/pi^2) * (2*dlw(alpha) + ddlw(alpha));
        diff = - ((2*alpha + adj_eps)/pi^2)*2*dlw(alpha) - ((2*alpha+1)/pi^2)*ddlw(alpha)*adj_eps;
        alpha = alpha - func/diff;
        
        crtn = abs(func/diff);
    end    
    
end
