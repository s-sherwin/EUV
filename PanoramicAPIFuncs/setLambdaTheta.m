function lambdaTheta = setLambdaTheta(lambdaTheta_full,ii,materials,nk,lambdaStr,thetaStr)
    %% Update the wavelength range and n/k
%     jj = [ii; ii];
    lambdaTheta = lambdaTheta_full(ii,:);
    setVariableValues(thetaStr,lambdaTheta(:,2));
    setVariableValues(lambdaStr,lambdaTheta(:,1));
    
    [a,b,c] = unique(lambdaTheta_full(:,1),'stable');
    kk = c(ii);
    %% Materials
    for iMat = 1:length(materials)
        n = real(nk(iMat,kk));
        setVariableValues(['n' materials{iMat}],n)
        k = -imag(nk(iMat,kk));
        setVariableValues(['k' materials{iMat}],k)
    end
end