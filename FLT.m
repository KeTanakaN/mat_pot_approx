function f = FLT(z) 
    if(z==0) 
        f = -pi^2/4;
    else
        f = -pi.*tanh(pi*z/4)./z;
    end
end
