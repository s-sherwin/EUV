function out = mkMLAbsEtch(n_cap,n_per,n_etch,thick,amu,rho_nom,rough,f0f1_elements,composition,lambdaTheta)
    out = [];
    N = size(thick,2);
    f = @(x) reshape(x,[],N);
    
    %% Find which layers are shared between etched and unetched ("mirror")
    if n_etch > 0
        inds_etch = 1:n_etch + 1;
    else
        inds_etch = [];
    end
    inds_detch = n_cap + n_per + 1;
    inds_contam = inds_detch(end)+1:size(thick,1);
    inds_shared = setdiff(1:size(thick,1),[inds_etch,inds_detch,inds_contam]);
    
    %% Mirror
    thick_mirror = f( thick(inds_shared,:) );
    composition_mirror = reshape( composition(:,inds_shared,:),size(composition,1),[],N );
    rough_mirror = f( rough(inds_shared,:) );
    n_capML = size(thick_mirror,1) - n_per;
    
    [layer_thick_ML,layer_nk_ML,layer_rough_ML] = mkMLStack(1,n_capML,thick_mirror,amu,rho_nom,rough_mirror,f0f1_elements,composition_mirror,lambdaTheta(:,1));
    out.mirror = [];
    out.mirror.thick = layer_thick_ML;
    out.mirror.nk = layer_nk_ML;
    out.mirror.rough = layer_rough_ML;
    out.mirror.rough(end,:) = 0;
    %% Absorber (un-etched)
    thick_abs = f( thick(inds_etch,:) );
    composition_abs = reshape( composition(:,inds_etch,:),size(composition,1),[],N );
    rough_abs = f( rough(inds_etch,:) );
    n_capAbs = size(thick_abs,1);
    
    [layer_thick_abs,layer_nk_abs,layer_rough_abs] = mkMLStack(0,n_capAbs,thick_abs,amu,rho_nom,rough_abs,f0f1_elements,composition_abs,lambdaTheta(:,1));
    out.absorber = [];
    out.absorber.thick = layer_thick_abs;
    out.absorber.nk = layer_nk_abs;
    out.absorber.rough = layer_rough_abs;
    out.absorber.rough(end,:) = 0;
    %% Etched absorber (+ contamination) 
    %% Top to bottom: Vacuum, contamination, remaining absorber, remaining cap
    detch = thick(inds_detch,:);
    % Bottom: Remaining layers from absorber
    thick_etch = thick_abs(n_etch:end,:);
    thick_etch(1,:) = detch;
    rough_etch = rough_abs(n_etch:end,:);
    rough_etch(1,:) = rough(inds_detch);
    composition_etch = composition_abs(:,end-size(thick_etch,1)+1:end,:);
    % On top of that: Carbon contamination
    thick_contam = thick(inds_contam,:);
    thick_etch = [thick_contam;thick_etch];
    composition_etch = [composition(:,inds_contam,:),composition_etch];
    rough_etch = [rough(inds_contam,:); rough_etch];
    % On top of that: Vacuum to match absorber reference plane
    delta_t = sum(thick_abs,1) - sum(thick_etch,1);
    thick_etch = [delta_t;thick_etch];
    rough_etch = [zeros(1,N); rough_etch];
    composition_etch = [0*composition_etch(:,1,:),composition_etch];
    
    if any(thick_etch < 0)
        is_neg = find(thick_etch < 0,1);
        t = thick_etch(is_neg,:);
        thick_etch(is_neg,:) = 0;
        thick_etch(is_neg+1,:) = thick_etch(is_neg+1,:) + t;
    end
    
    thick_etch = f(thick_etch);
    composition_etch = reshape( composition_etch,size(composition,1),[],N );
    rough_etch = f(rough_etch);
    
    [layer_thick_etch,layer_nk_etch,layer_rough_etch] = mkMLStack(0,size(thick_etch,1),thick_etch,amu,rho_nom,rough_etch,f0f1_elements,composition_etch,lambdaTheta(:,1));
    out.etch = [];
    out.etch.thick = layer_thick_etch;
    out.etch.nk = layer_nk_etch;
    out.etch.rough = layer_rough_etch;
    out.etch.rough(end,:) = 0;
end