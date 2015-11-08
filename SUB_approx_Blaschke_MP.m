function fapp = SUB_approx_Blaschke_MP(x, sample, func , LogInvWeight)
    a = sample;
    f = func;
    liw = LogInvWeight;

    %% Setting an x-grid and loading the sampling points
    lx = length(x);
    la = length(a);

    %% Approximation of a function
    LIW = mp(zeros(1,lx)); % MULTI
    for l=1:lx
        LIW(l) = liw(x(l));
    end

    Bprod_x = mp(ones(la,lx)); % MULTI
    Bprod_a = mp(ones(1,la)); % MULTI
    for j=1:la
        for n=1:la
            if(n~=j)
                Bprod_x(j,:) = Bprod_x(j,:) .* tanh(x-a(n));
                Bprod_a(j) = Bprod_a(j) * tanh(a(j)-a(n));
            end
        end
        Bprod_x(j,:) = Bprod_x(j,:) .* exp(-LIW) ./ (cosh(x-a(j)).^2);
        Bprod_a(j) = Bprod_a(j) * exp(-liw(a(j)));
    end
    fapp = (f(a) ./ Bprod_a) * Bprod_x;
end
