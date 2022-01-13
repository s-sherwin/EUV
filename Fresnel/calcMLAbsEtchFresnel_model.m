function outputs = calcMLAbsEtchFresnel_model(model,illum_tbl)
    materials = model.materials;
    ind_struct = model.ind_struct;
    x = model.x;
    n_etch = model.n_etch;
    n_per = model.n_per;
    N_ML = model.N_ML;
    n_cap = model.n_cap;
    pol = model.pol;
    n_contam = model.n_contam;
    %%
    composition = x(ind_struct.composition);
    composition = reshape(composition,length(materials),[]);
    lambdaTheta1 = illum_tbl(:,1:2);
    [f0f1_elements,rho_nom,amu] = getScatFac(materials,lambdaTheta1(:,1));

    fFres = @(x) calcMLAbsEtchFresnel(N_ML,n_cap,n_per,n_contam,n_etch,x(ind_struct.thick,:),amu,rho_nom,x(ind_struct.rough,:),f0f1_elements,reshape(x(ind_struct.composition,:),[size(composition),size(x,2)]),lambdaTheta1,pol,1,x(ind_struct.substrate_roughness,:));
    [r_MLabs,r_MLetch,r_ML,r_abs,t_abs,r_etch,t_etch,M_ML,M_abs,M_etch,M_MLabs,M_MLetch] = fFres(x);
    %%
    outputs = [];
    outputs.r = cat(2,r_MLabs,r_MLetch,r_ML,r_abs,r_etch);
    outputs.t = cat(2,t_abs,t_etch);
    outputs.M = cat(4,M_ML,M_abs,M_etch,M_MLabs,M_MLetch);
end