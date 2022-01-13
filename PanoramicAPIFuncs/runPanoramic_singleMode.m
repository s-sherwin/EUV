function outStruct = runPanoramic_singleMode(simParams,illumParams,materialParams,geometryParams,blankStr,runType,nPlanes,saveLoc)
    if nargin < 7
        nPlanes = 1;
    end
    simFields = fieldnames(simParams);
    illumFields = fieldnames(illumParams);
    if nargin < 8
        saveLoc = '';
    end
    %% Loop through and run simulations
    if isfield(simParams,'p_y') % 3D XYZ sim
        sim_dim = 3;
        px = simParams.mag_x*simParams.p_x;
        py = simParams.mag_y*simParams.p_y;
        N = floor(2*[px py]/illumParams.lambda);
        orders_X = -N(1):N(1);
        orders_Y = -N(2):N(2);
        scatteredOrders = zeros(length(orders_Y),length(orders_X),2,2,nPlanes);
        if isfield(simParams,'dx')
            Nx = simParams.p_x/simParams.dx;
            Ny = simParams.p_y/simParams.dx;
%         elseif simParams.p_x ~= simParams.p_y
%             Nx = simParams.Nx;
%             Ny = round(simParams.Nx*p_y/p_x); % Note: must round to have an integer number of pixels; p_y should usually be equal to p_x or double for anamorphic magnification
        else
            Nx = simParams.Nx;
            Ny = round(Nx*py/px); % Note: must round to have an integer number of pixels; p_y should usually be equal to p_x or double for anamorphic magnification
        end
        eField = zeros(Ny,Nx,2,3,nPlanes);
    else % 2D XZ sim
        sim_dim = 2;
        N = floor(2*simParams.p*simParams.mag/illumParams.lambda);
        orders = -N:N;
        scatteredOrders = zeros(length(orders),2,2,nPlanes);
        eField = zeros(simParams.Nx,2,3,nPlanes);
    end
    %% Load simulation with blank stack
    loadSetup(blankStr);
    %% Set simulation parameters
    for j = 1:length(simFields)
        setVariableValues(simFields{j},simParams.(simFields{j}));
    end
    %% Set illumination parameters for current settings
    for j = 1:length(illumFields)
        setVariableValues(illumFields{j},illumParams.(illumFields{j}));
    end
    %% Set material properties for current settings
    for j = 1:length(materialParams.materials)
        str = ['n' materialParams.materials{j}];
        if isfield(materialParams,str)
            setVariableValues(str,materialParams.(str))
        else
            setVariableValues(str,real(materialParams.nk(1,j)));
        end
        str = ['k' materialParams.materials{j}];
        if isfield(materialParams,str)
            setVariableValues(str,materialParams.(str))
        else
            setVariableValues(str,-imag(materialParams.nk(1,j)));
        end
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
                if t > 0
                    coordString = coords2String([z0 z0 + t]);
                    z0 = z0 + t;
                    callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
                end
            end
        end
        Cap = geometryParams.stack.Cap;
        for k = length(Cap.t):-1:1
            material = Cap.material{k};
            t = Cap.t(k);
            if t > 0
                coordString = coords2String([z0 z0 + t]);
                z0 = z0 + t;
                callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
            end
        end
    end
    %% Geometry: background material (if not vacuum)
    if isfield(geometryParams,'backgroundMatStr')
        t = geometryParams.t;
        material = geometryParams.backgroundMatStr;
        coordString = coords2String([z0 z0 + t]);
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
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
    elseif strcmp(geometryParams.type,'box2D')
        coordsX = geometryParams.coordsX;
        coordsY = geometryParams.coordsY;
        Dx = geometryParams.Dx;
        Dy = geometryParams.Dy;
        t = geometryParams.t;
        matStr = geometryParams.matStr;

        coordsX_tmp = (coordsX - mean(coordsX))*Dx + mean(coordsX);
        coordsY_tmp = (coordsY - mean(coordsY))*Dy + mean(coordsY);
        coords = [coordsX_tmp coordsY_tmp z0 z0+t];
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
    elseif strcmp(geometryParams.type,'symTrapStack')
        matStr = geometryParams.matStr;
        coordsX = geometryParams.coordsX;
        coordsY = geometryParams.coordsY;
        
        t = geometryParams.t;
        x = geometryParams.x*coordsX(2); 
        z = geometryParams.z;
        if length(x) == length(z) + 2 % 2 additional x coords than z coords
            z = [0; vec(z); t];
        end
        
        for i = 1:length(x)-1
            X = coordsX(2)/2 + [-x(i) -x(i+1) x(i+1) x(i)];
%             X = x0*[1 - wBottom, 1 - wTop, 1 + wTop, 1 + wBottom];
            dz = z(i+1) - z(i);
            Z = [z0, z0 + dz, z0 + dz, z0];
            z0 = z0 + dz;
            coords = [coordsY, X(1), Z(1), X(2), Z(2), X(3), Z(3), X(4), Z(4)];%, Z(1), X(1)];
            coordString = coords2String(coords);
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','yPolygon',matStr,coordString);
        end
    elseif strcmp(geometryParams.type,'sphere')
        t = geometryParams.t;
        matStr = geometryParams.matStr;
        x0 = geometryParams.x0;
        y0 = geometryParams.y0;
        z0 = geometryParams.z0;
        R = geometryParams.R;
        coords = [x0 - R, x0 + R, y0 - R, y0 + R, z0 - R, z0 + R];
        coordString = coords2String(coords);
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','yPolygon',matStr,coordString);
    end
    z0 = z0 + t;
    %%
    if isfield(geometryParams,'box2')
        coordsX = geometryParams.coordsX;
        coordsY = geometryParams.coordsY;
        D = geometryParams.D;
        t = geometryParams.box2.t;
        dz = geometryParams.box2.dz;
        z0 = z0 + dz;
        
        if isfield(geometryParams,'backgroundMatStr')
            material = geometryParams.backgroundMatStr;
            coordString = coords2String([z0 z0 + t]);
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','layer',material,coordString);
        end
        
        if isfield(geometryParams.box2,'matStr')
            matStr = geometryParams.box2.matStr;
        else
            matStr = geometryParams.matStr;
        end
        coordsX_tmp = (coordsX - mean(coordsX))*D + mean(coordsX);
        coords = [coordsX_tmp coordsY z0 z0+t];
        coordString = coords2String(coords);
        material = matStr;
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
    end
    %%
    if ~isempty(saveLoc)
        [filepath] = fileparts(saveLoc);
        if exist(filepath,'dir')
            %% Save setup
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','saveSetup',saveLoc)
%             callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','saveSetup','C:\Users\stuar\Documents\Panoramic\Scatterometry\test.sim')
        end
    end

    %% Run simulation
    if strcmp(runType,'FDTD')
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingTEMPESTpr2',1);
    elseif strcmp(runType,'RCWA')
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingRCWA',1);
    elseif strcmp(runType,'TMM')
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingStackKirchhoffThick',1);
    end

    %% Imaging calculation
    try
        javaMethod('setObject','panoramictech.v700.OpenAPI.CAPIClient',outh)
        outh2 = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateImage',1);
        [table,header] = getEMOutputTable(outh2,'Resist Image');
        if sim_dim == 2
            if length(header) == 3
                x_AI = unique(table(:,1));
                z_AI = unique(table(:,2));
                AI = reshape(table(:,3),[],length(x_AI));
            end
        else

        end
    catch
        x_AI = [];
        z_AI = [];
        AI = [];
    end
    %% Collect scattered orders
    [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
    [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);

    if sim_dim == 2
        %% Organize simulation outputs
        tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        cgrp = unique(tableTE(:,strcmp(header,'cgrp')));
        %% Pull out the orders 
        for iCGRP = 1:length(cgrp)
            col_ind = 3;
            rows = tableTE(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            if length(header) > 3
                rows = and(rows,tableTE(:,strcmp(header,'n')) == 0);
                col_ind = strcmp(header,'Data');
            end
            %% TE Output
            scatteredOrders(:,iCGRP,1,:) = reshape(tableTE(rows,col_ind),[],1,1,nPlanes);
            %% TM Output
            scatteredOrders(:,iCGRP,2,:) = reshape(tableTM(rows,col_ind),[],1,1,nPlanes);
        end

        %% Collect E-fields
        [tableEx,~] = getEMOutputTable(outh,'Ex',true);
        [tableEy,~] = getEMOutputTable(outh,'Ey',true);
        [tableEz,header] = getEMOutputTable(outh,'Ez',true);

        xx = unique(tableEx(:,1));

        for iCGRP = 1:length(cgrp)
            rows = tableEx(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            eField(:,iCGRP,1,:) = reshape(tableEx(rows,3),[],1,1,nPlanes);
            eField(:,iCGRP,2,:) = reshape(tableEy(rows,3),[],1,1,nPlanes);
            eField(:,iCGRP,3,:) = reshape(tableEz(rows,3),[],1,1,nPlanes);
        end
    else
        %% Organize simulation outputs
        tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTE(isnan(tableTE(:,strcmp(header,'m'))),strcmp(header,'m')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'m'))),strcmp(header,'m')) = 0;
        cgrp = unique(tableTE(:,strcmp(header,'cgrp')));
        %%
        nx = (length(unique(round(tableTE(:,strcmp(header,'m')),1)))-1)/2;
        ny = (length(unique(round(tableTE(:,strcmp(header,'n')),1)))-1)/2;
        ii_x = ismember(-nx:nx,orders_X);
        ii_y = ismember(-ny:ny,orders_Y);
        if isempty(ii_y)
            ii_y = 1;
        end
        if isempty(ii_x)
            ii_x = 1;
        end
        %% Pull out the orders 
        for iCGRP = 1:length(cgrp)
            %% TE Output
            rows = tableTE(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            tmp = reshape(tableTE(rows,strcmp(header,'Data')),[],length(ii_x),1,1,nPlanes);
            scatteredOrders(:,:,iCGRP,1,:) = tmp(ii_y,ii_x,:,:,:);
%             nx = length(unique(round(tableTE(:,strcmp(header,'m')),1)));
%             scatteredOrders(:,:,iCGRP,1,:) = reshape(tableTE(rows,strcmp(header,'Data')),[],nx,1,1,nPlanes);
            %% TM Output
            tmp = reshape(tableTM(rows,strcmp(header,'Data')),[],length(ii_x),1,1,nPlanes);
            scatteredOrders(:,:,iCGRP,2,:) = tmp(ii_y,ii_x,:,:,:);
%             scatteredOrders(:,:,iCGRP,2,:) = reshape(tableTM(rows,strcmp(header,'Data')),[],length(orders_X),1,1,nPlanes);
        end

        %% Collect E-fields
        [tableEx,~] = getEMOutputTable(outh,'Ex',true);
        [tableEy,~] = getEMOutputTable(outh,'Ey',true);
        [tableEz,header] = getEMOutputTable(outh,'Ez',true);
        
        xx = unique(tableEx(:,1));
        yy = unique(tableEx(:,2));

        for iCGRP = 1:length(cgrp)
            rows = tableEx(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            eField(:,:,iCGRP,1,:) = reshape(tableEx(rows,4),[],length(xx),1,1,nPlanes);
            eField(:,:,iCGRP,2,:) = reshape(tableEy(rows,4),[],length(xx),1,1,nPlanes);
            eField(:,:,iCGRP,3,:) = reshape(tableEz(rows,4),[],length(xx),1,1,nPlanes);
        end
    end
    %% Modulate negative orders by following phase filters:
    if sim_dim == 2
       %%  shift negative angles by pi
       deltaST = illumParams.lambda/(simParams.p*simParams.mag);
       sinT = sind(illumParams.theta) + orders*deltaST;
       jj = sinT < 0;
       %%  shift odd orders by pi
       scatteredOrders = scatteredOrders.*exp(1i*pi*jj').*exp(1i*pi*(mod(orders',2)==1));
%         else % Original code for on-axis
%             %%  shift negative orders by pi
%             %%  shift odd orders by pi
%             scatteredOrders = scatteredOrders.*exp(1i*pi*(orders < 0)).*exp(1i*pi*(mod(orders,2)==1));
%         end
    else
        if false
           %%  shift negative angles by pi
           deltaST_x = illumParams.lambda/simParams.p_x;
           deltaST_y = illumParams.lambda/simParams.p_y;
           sinT_x = sind(illumParams.theta)*cosd(illumParams.phi) + orders_X*deltaST_x;
           jj = sinT_x < 0;
           sinT_y = sind(illumParams.theta)*sind(illumParams.phi) + orders_Y*deltaST_y;
           kk = sinT_y < 0;
           %%  shift odd orders by pi
           scatteredOrders = scatteredOrders.*exp(1i*pi*jj).*exp(1i*pi*kk').*exp(1i*pi*(mod(orders_X + orders_Y',2)==1));
        end
    end
    
    %% Prepare outputs
    outStruct.scatteredOrders = scatteredOrders;
    outStruct.eField = eField;
    if exist('AI','var')
        outStruct.AI = AI;
    end
    if sim_dim == 2
        outStruct.orders = orders;
        outStruct.xx = xx;
        if exist('x_AI','var')
            outStruct.xAI = x_AI;
        end
        if exist('z_AI','var')
            outStruct.zAI = z_AI;
        end
    else
        outStruct.orders_X = orders_X;
        outStruct.orders_Y = orders_Y;
        outStruct.xx = xx;
        outStruct.yy = yy;
    end
    
    %% Delete results for memory management
    destroyAll();
end