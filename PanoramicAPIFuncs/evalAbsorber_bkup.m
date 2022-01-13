function [ScatteredOrders,outh] = evalAbsorber(blankStr,Nx,Nz,stack,materials,nk,coordsX,coordsY,...
    lambdaTheta_full,ii,lambdaStr,thetaStr,tAbsorber,wTopAbsorber,wBottomAbsorber,matAbsorber,runStr,orders)
    %% Load simulation with blank stack
    loadSetup(blankStr);

    %% Update the grid size
    dz = stack.Periodic.dspace/Nz;
    setVariableValues('dz',dz);
    setVariableValues('Nx',Nx);
    
    %% Add a multilayer stack
    z0 = createUnpatternedStack_API(stack,coordsX,coordsY);

    %% Set material properties and nk
    lambdaTheta = setLambdaTheta(lambdaTheta_full,ii,materials,nk,lambdaStr,thetaStr);

    %% Add absorber 
    if wTopAbsorber == wBottomAbsorber
        %% Add recctangular absorber block
        coordsXAbsorber = mean(coordsX) + (coordsX - mean(coordsX))*wTopAbsorber;
        coords = [coordsXAbsorber coordsY z0 z0+tAbsorber];
        coordString = coords2String(coords);
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',matAbsorber,coordString);
    else
        %% Add trapezoidal absorber block
        x0 = mean(coordsX);
        X = x0*[1 - wBottomAbsorber, 1 - wTopAbsorber, 1 + wTopAbsorber, 1 + wBottomAbsorber];
%         X = [x0 - wBottomAbsorber/2, x0 - wTopAbsorber/2, x0 + wTopAbsorber/2, x0 + wBottomAbsorber/2];
        Z = [z0, z0 + tAbsorber, z0 + tAbsorber, z0];
%         coords = [coordsY, Z(1), X(1), Z(2), X(2), Z(3), X(3), Z(4), X(4), Z(1), X(1)];
        coords = [coordsY, X(1), Z(1), X(2), Z(2), X(3), Z(3), X(4), Z(4)];%, Z(1), X(1)];
        coordString = coords2String(coords);
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','yPolygon',matAbsorber,coordString);
    end
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
    [~,~,cols] = intersect({lambdaStr,'cgrp','m','n'},header,'stable');
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
    rows = ismember(tableTM(:,cols),outCols,'rows');
    R_TM = tableTM(rows,dataCol);

%     for iCGRP = 1:length(cgrp)
%         for iPol = 1:2
%             if iPol == 1 % TE
%                 
%             else % TM
%                 
%             end
%         end
%         R_hyp_pol(iCGRP,iPol) = tmp;
%     end
    R_hyp_pol(:,:,1,1) = reshape(abs(R_TE(1:end/2)).^2,length(m),[])';
    R_hyp_pol(:,:,2,1) = reshape(abs(R_TE(end/2+(1:end/2))).^2,length(m),[])';
    R_hyp_pol(:,:,1,2) = reshape(abs(R_TM(1:end/2)).^2,length(m),[])';
    R_hyp_pol(:,:,2,2) = reshape(abs(R_TM(end/2+(1:end/2))).^2,length(m),[])';
    
    s = size(R_hyp_pol);
    ScatteredOrders = sum(reshape(R_hyp_pol,s(1),s(2),[]),3);
end