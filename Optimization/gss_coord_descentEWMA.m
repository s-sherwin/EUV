function [x_fit,fx] = gss_coord_descentEWMA(x0,f,LB,UB,ncycles,niter,fac,alpha,step0)
    if nargin < 9
        step0 = 1;
    end
    x_fit = x0;
    L = LB;
    U = UB;
    D = U - L;
    D0 = D*step0;
    for iC = 1:ncycles
        %% Coordinate descent
        for i = 1:length(x0)
            ftmp = @(v) f(f_ii(v,x_fit,(1:length(x0))==i));
            %% Golden-Section Search
            dx = D0(i)*fac^(iC-1);
            if abs(dx) < 1e-14
                continue
            end
            L1 = max(L(i),x_fit(i)-dx);
            U1 = min(U(i),x_fit(i)+dx);
            [v,c] = gss(L1,U1,@(x) ftmp(x),niter);
            %% Update with EWMA
            x_fit(i) = v*alpha + x_fit(i)*(1 - alpha);
        end
        %% Update bounds
        D = U - L;
        D = max(D*fac,D0*fac^(iC-1));
        L = x_fit - D/2;
        L = max(L,LB);
        U = x_fit + D/2;
        U = min(U,UB);
    end
    if nargout > 1
        fx = f(x_fit);
    end
end