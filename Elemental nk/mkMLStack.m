function [layer_thick,layer_nk,layer_rough] = mkMLStack(N_ML,n_cap,thick,amu,rho_nom,rough,f0f1,composition,lambda)
    if isempty(thick)
        layer_thick = [];
        layer_nk = [];
        layer_rough = [];
        return
    end
    ii_cap = 1:n_cap;
    ii_per = n_cap + 1:size(thick,1);
    layer_thick = [thick(ii_cap,:); repmat(thick(ii_per,:),N_ML,1)];
    layer_rough = [rough(ii_cap,:); repmat(rough(ii_per,:),N_ML,1); rough(end,:)];
    
    re = 2.8179e-15; % classical electron radius [m]
    Na = 6.0221409e23; % Avogodro's number
    na = rho_nom./amu * Na; % atoms/cm^3
    const = 1e-12 * lambda.^2 .* re / (2*pi);
    
    composition = reshape(composition,[],size(thick,1),size(thick,2));
    
    nk = 1 - const .*  f0f1*reshape(na.*composition,length(na),[]);
    nk = reshape(nk,length(lambda),size(thick,1),size(thick,2));
    layer_nk = [nk(:,ii_cap,:) repmat(nk(:,ii_per,:),1,N_ML,1)];
end