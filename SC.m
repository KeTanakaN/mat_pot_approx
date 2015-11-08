function f = SC(alpha, v) 
    if(v==0) 
        f = alpha + 1/2;
    else
        f = 2*(2*sin(alpha*v) + v.*cos(alpha*v))./(v.*(4+v.^2));
    end
end
