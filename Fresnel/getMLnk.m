function [layer_nk] = getMLnk(N_ML,n_cap,n_per,amu,rho_nom,f0f1,composition,lambda)
    n = n_cap + n_per;
    ii_cap = 1:n_cap;
    ii_per = n_cap + (1:n_per);
    
    composition = reshape(composition,length(rho_nom),n);
    
    re = 2.8179e-15; % classical electron radius [m]
    Na = 6.0221409e23; % Avogodro's number
    na = rho_nom./amu * Na; % atoms/cm^3
    const = 1e-12 * lambda.^2 .* re / (2*pi);
    
%     composition = reshape(composition,[],size(thick,1),size(thick,2));
    
    nk = 1 - const .*  f0f1*reshape(na.*composition,length(na),[]);
    nk = reshape(nk,length(lambda),n);
    layer_nk = [nk(:,ii_cap) repmat(nk(:,ii_per),1,N_ML,1)];
end