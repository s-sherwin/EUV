function out_mat = calcMLAbsEtchFresnel_wrapper(X,fFres)
    [r_MLabs,r_MLetch,r_ML,r_abs,t_abs,r_etch,t_etch,M_ML,M_abs,M_MLabs,M_MLetch] = fFres(X);
    out_mat = cat(3,r_MLabs,r_MLetch,r_ML,r_abs,t_abs,r_etch,t_etch);
end