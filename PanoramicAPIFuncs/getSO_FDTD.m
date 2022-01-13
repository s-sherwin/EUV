function scatteredOrders = getSO_FDTD(simParams,illumParams,materialParams,geometryParams,blankStr,orders,nPlanes)
    if nargin < 7
        nPlanes = 1;
    end
    simFields = fieldnames(simParams);
    illumFields = fieldnames(illumParams);
    materialFields = fieldnames(materialParams);
    geometryFields = fieldnames(geometryParams);
    %% Loop through and run simulations
    scatteredOrders = zeros(length(illumParams.(illumFields{1})),length(orders),2,2,nPlanes);
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
                for i = length(Periodic.t):-1:1
                    material = Periodic.material{i};
                    t = Periodic.t(i)*Periodic.dspace;
                    coordString = coords2String([z0 z0 + t]);
                    z0 = z0 + t;
                    callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
                end
            end
            Cap = geometryParams.stack.Cap;
            for i = length(Cap.t):-1:1
                material = Cap.material{i};
                t = Cap.t(i);
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
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingTEMPESTpr2',1);

        %% Collect outputs
        [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
        [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);
        
        %% Organize simulation outputs
        m = 0.5 + orders;
        tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        cgrp = unique(tableTE(:,strcmp(header,'cgrp')));

        [~,~,cols] = intersect({'cgrp','m'},header,'stable');
        dataCol = strcmp(header,'Data');

        %% Pull out the 0 order for each angle
        cgrp = [0 1];
        for iCGRP = 1:2
            for iO = 1:length(orders)
                %% TE Output
                rows = ismember(tableTE(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(i,iO,iCGRP,1,:) = tableTE(rows,dataCol);

                %% TM Output
                rows = ismember(tableTM(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(i,iO,iCGRP,2,:) = tableTM(rows,dataCol);
            end
        end
    end
    
    %% Modulate negative orders by following phase filters:
    %%  shift negative orders by pi
    %%  shift odd orders by pi
    scatteredOrders = scatteredOrders.*exp(1i*pi*(orders < 0)).*exp(1i*pi*(mod(orders,2)==1));
end