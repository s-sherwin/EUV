function [scatteredOrders,eField,I] = runPanoramic(simParams,illumParams,materialParams,geometryParams,blankStr,orders,runType,nPlanes)
    if nargin < 8
        nPlanes = 1;
    end
    simFields = fieldnames(simParams);
    illumFields = fieldnames(illumParams);
    materialFields = fieldnames(materialParams);
    geometryFields = fieldnames(geometryParams);
    %% Loop through and run simulations
    scatteredOrders = zeros(length(illumParams.(illumFields{1})),length(orders),2,2,nPlanes);
    eField = zeros(length(illumParams.(illumFields{1})),simParams.Nx,2,3,nPlanes);
    I = zeros(length(illumParams.(illumFields{1})),simParams.Nx,2,nPlanes);
    for i = 1:length(illumParams.(illumFields{1}))
        %% Load simulation with blank stack
        loadSetup(blankStr);
        %% Set simulation parameters
        for j = 1:length(simFields)
            setVariableValues(simFields{j},simParams.(simFields{j}));
        end
        %% Set illumination parameters for current settings
        for j = 1:length(illumFields)
            setVariableValues(illumFields{j},illumParams.(illumFields{j})(i));
        end
        %% Set material properties for current settings
        for j = 1:length(materialParams.materials)
            setVariableValues(['n' materialParams.materials{j}],real(materialParams.nk(i,j)));
            setVariableValues(['k' materialParams.materials{j}],-imag(materialParams.nk(i,j)));
        end
        %% Set geometry
        z0 = geometryParams.z0;
        %% Geometry: multilayer stack
        if isfield(geometryParams,'stack')
            %% Add layers for multiayer stack
            Periodic = geometryParams.stack.Periodic;
            for ind = 1:Periodic.N
                for k = length(Periodic.t):-1:1
                    material = Periodic.material{k};
                    t = Periodic.t(k)*Periodic.dspace;
                    coordString = coords2String([z0 z0 + t]);
                    z0 = z0 + t;
                    callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
                end
            end
            Cap = geometryParams.stack.Cap;
            for k = length(Cap.t):-1:1
                material = Cap.material{k};
                t = Cap.t(k);
                coordString = coords2String([z0 z0 + t]);
                z0 = z0 + t;
                callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
            end
        end
        %% Geometry: absorber
        if strcmp(geometryParams.type,'box')
            coordsX = geometryParams.coordsX;
            coordsY = geometryParams.coordsY;
            D = geometryParams.D;
            t = geometryParams.t;
            matStr = geometryParams.matStr;
            
            coordsX_tmp = (coordsX - mean(coordsX))*D + mean(coordsX);
            coords = [coordsX_tmp coordsY z0 z0+t];
            coordString = coords2String(coords);
            material = matStr;
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
        elseif strcmp(geometryParams.type,'symTrap')
            coordsX = geometryParams.coordsX;
            coordsY = geometryParams.coordsY;
            D = geometryParams.D;
            t = geometryParams.t;
            deltaX = geometryParams.deltaX;
            matStr = geometryParams.matStr;
            
            x0 = mean(coordsX);
            wBottom = D + deltaX/2;
            wTop = D - deltaX/2;
            X = x0*[1 - wBottom, 1 - wTop, 1 + wTop, 1 + wBottom];
            Z = [z0, z0 + t, z0 + t, z0];
            coords = [coordsY, X(1), Z(1), X(2), Z(2), X(3), Z(3), X(4), Z(4)];%, Z(1), X(1)];
            coordString = coords2String(coords);
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','yPolygon',matStr,coordString);
        end
        
        %%
        if false
            %% Save setup
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','saveSetup','C:\Users\stuar\Documents\Panoramic\Scatterometry\test.sim')
        end
        
        %% Run simulation
        if strcmp(runType,'FDTD')
            outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingTEMPESTpr2',1);
        elseif strcmp(runType,'RCWA')
            outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingRCWA',1);
        end
        
        %% Collect scattered orders
        [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
        [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);
        
        %% Organize simulation outputs
        m = 0.5 + orders;
        tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        cgrp = unique(tableTE(:,strcmp(header,'cgrp')));

        [~,~,cols] = intersect({'cgrp','m'},header,'stable');
        dataCol = strcmp(header,'Data');

        %% Pull out the orders 
%         cgrp = [0 1];
        for iCGRP = 1:length(cgrp)
            for iO = 1:length(orders)
                %% TE Output
                rows = ismember(tableTE(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(i,iO,iCGRP,1,:) = tableTE(rows,dataCol);

                %% TM Output
                rows = ismember(tableTM(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(i,iO,iCGRP,2,:) = tableTM(rows,dataCol);
            end
        end
        
        %% Collect E-fields
        [tableEx,~] = getEMOutputTable(outh,'Ex',true);
        [tableEy,~] = getEMOutputTable(outh,'Ey',true);
        [tableEz,header] = getEMOutputTable(outh,'Ez',true);
        
        for iCGRP = 1:length(cgrp)
            rows = tableEx(:,2) == cgrp(iCGRP);
            eField(i,:,iCGRP,1,:) = reshape(tableEx(rows,3),1,[],1,1,length(nPlanes));
            eField(i,:,iCGRP,2,:) = reshape(tableEy(rows,3),1,[],1,1,length(nPlanes));
            eField(i,:,iCGRP,3,:) = reshape(tableEz(rows,3),1,[],1,1,length(nPlanes));
        end
        %% Collect intensities
        [tableI,header] = getEMOutputTable(outh,'Electric Field Intensity',false);
        
        for iCGRP = 1:length(cgrp)
            rows = tableI(:,2) == cgrp(iCGRP);
            I(i,:,iCGRP,:) = reshape(tableI(rows,3),1,[],1,1,length(nPlanes));
        end
    end
    
    %% Modulate negative orders by following phase filters:
    if true
       for i = 1:length(illumParams.(illumFields{1}))
           %%  shift negative angles by pi
           deltaST = illumParams.lambda(i)/(simParams.p*simParams.mag);
           sinT = sind(illumParams.theta(i)) + orders*deltaST;
           jj = sinT < 0;
           %%  shift odd orders by pi
           scatteredOrders = scatteredOrders.*exp(1i*pi*jj).*exp(1i*pi*(mod(orders,2)==1));
       end
    else % Original code for on-axis
        %%  shift negative orders by pi
        %%  shift odd orders by pi
        scatteredOrders = scatteredOrders.*exp(1i*pi*(orders < 0)).*exp(1i*pi*(mod(orders,2)==1));
    end
end