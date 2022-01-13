function [ScatteredOrders,outh] = evalEtchedPSM(stack,coordsX,coordsY,...
    lambdaTheta_full,ii,lambdaStr,thetaStr,NEtch,wEtch,runStr,orders)
    %% Add a multilayer stack
    z0 = createUnpatternedStack_API(stack,coordsX,coordsY);

    %% Set material properties and nk
    lambdaTheta = setLambdaTheta(lambdaTheta_full,ii,stack,lambdaStr,thetaStr);

    %% Etch 
    etchDepth = sum(stack.Cap.t) + NEtch*stack.Periodic.dspace;
    coordsXEtch = mean(coordsX) + (coordsX - mean(coordsX))*wEtch;
    coords = [coordsXEtch coordsY z0-etchDepth z0+10];
    coordString = coords2String(coords);
    material = 'Vacuum';
    callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
    
    %% Run scattering simulation
    outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient',runStr,1); % Run scattering simulation
    
    %% Collect outputs
    [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
    [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);
    %% Organize simulation outputs
    m = 0.5 + orders;
    tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
    tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
    cgrp = unique(tableTE(:,strcmp(header,'cgrp')));
    R_hyp_pol = zeros(size(lambdaTheta,1),length(m),length(cgrp),2);
    [~,~,cols] = intersect({'wiz_wavelength','cgrp','m','n'},header,'stable');
    dataCol = strcmp(header,'Data');
    %%
    outCols = zeros(0,length(cols));
    for iLam = 1:size(lambdaTheta,1)
        for iOrder = 1:length(m)
            for iCGRP = 1:length(cgrp)
                %%
                outCols = [outCols; lambdaTheta(iLam,1) cgrp(iCGRP) m(iOrder) zeros(1,length(cols) - length(header) + 1)];
            end
        end
    end
    %%
    rows = ismember(tableTE(:,cols),outCols,'rows');
    R_TE = tableTE(rows,dataCol);
    R_TM = tableTM(rows,dataCol);

    R_hyp_pol(:,:,1,1) = reshape(abs(R_TE(1:end/2)).^2,length(m),[])';
    R_hyp_pol(:,:,2,1) = reshape(abs(R_TE(end/2+(1:end/2))).^2,length(m),[])';
    R_hyp_pol(:,:,1,2) = reshape(abs(R_TM(1:end/2)).^2,length(m),[])';
    R_hyp_pol(:,:,2,2) = reshape(abs(R_TM(end/2+(1:end/2))).^2,length(m),[])';
    
    s = size(R_hyp_pol);
    ScatteredOrders = sum(reshape(R_hyp_pol,s(1),s(2),[]),3);
end