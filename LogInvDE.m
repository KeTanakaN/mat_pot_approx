function y = LogInvDE(x)
% In the case |x|>=6, 
% this function returns an approximate value of log(w_{DE}^{-1}(x)).
    if(abs(x)<6)
        y1 = log(1+exp(pi.*sinh(x)));
        y2 = log(1+exp(-pi.*sinh(x)));
    elseif(x>=6)
        y1 = pi.*sinh(x);
        y2 = 0;
    elseif(x<=-6)
        y1 = 0;
        y2 = -pi.*sinh(x);
    end
    y = (y1 + y2)/2.0 - log(2);
end
