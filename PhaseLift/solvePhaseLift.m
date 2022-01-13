function [E,X,cc,err_true] = solvePhaseLift(X0,b,B,Bpinv,Weight,settings)
    if nargin < 5
        Weight = ones(size(b));
    end
    if nargin < 6
        settings = [];
    end
    if isfield(settings,'alpha_lr')
        alpha_lr = settings.alpha_lr;
    else
        alpha_lr = 1e-2;
    end
    if isfield(settings,'step_mag')
        step_mag = settings.step_mag;
    else
        step_mag = 1e-8;
    end
    if isfield(settings,'step_decay')
        step_decay = settings.step_decay;
    else
        step_decay = 1 - 1e-2;
    end
    if isfield(settings,'N_GD')
        N_GD = settings.N_GD;
    else
        N_GD = 5;
    end
    if isfield(settings,'N_iter')
        N_iter = settings.N_iter;
    else
        N_iter = 40;
    end
    if isfield(settings,'Xtrue')
        Xtrue = settings.Xtrue;
    else
        Xtrue = 0*X0;
    end
    %%
    X = X0;
    cc = zeros(N_iter,N_GD);
    err_true = zeros(N_iter,N_GD);
    for i = 1:N_iter
        for j = 1:N_GD
            y = B*vec(X);
            y = mean_pad(y,Weight);
            res = y - b;
            res = res.*Weight;

            cc(i,j) = rms(res)/rms(b);

            step = Bpinv*res;
            step = step/norm(step)*step_mag*step_decay^(i-1);
            step = reshape(step,size(X));
            X = X - step;
            err_true(i,j) = norm(X - Xtrue,'fro');
        end
        %% Low-rank
        [U,S,V] = svd(X,'econ');
        s = diag(S);
        s = softThresh(s,alpha_lr*s(1));
        S = diag(s);
        X = U*S*U';
        %% Extract rank 1 component (recovered field)
        E = U(:,1);
        i0 = floor(length(E)/2)+1;
        E = E/exp(1i*angle(E(i0)));
        E = E/norm(E);
        %% Trace
        X = X/trace(X);
    end
end