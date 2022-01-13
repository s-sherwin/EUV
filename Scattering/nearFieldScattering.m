function E = nearFieldScattering(settings)
    %% Define the thin-film model (absorber, multilayer, etch) - required field
    model = settings.model; % Thin-film model
    %% Which scattering model to use - optional field
    if isfield(settings,'model_str')
        model_str = settings.model_str;
    else
        model_str = 'DblSc'; % Default scattering model
    end
    %% Define illumination settings - required field
    illum_tbl = settings.illum_tbl; % [lambda, theta, phi]; [nm, degress, degrees]
    %% Define pitch in x and y - required fields
    px = settings.px;
    py = settings.py;
    %% Define diffraction sampling in x and y - required fields
    orders_x = settings.orders_x;
    orders_y = settings.orders_y;
    %% Maximum angle - optional field
    if isfield(settings,'theta_max')
        theta_max = settings.theta_max;
    else
        theta_max = 80; % Default maximum angle [degrees]
    end
    %% Low-pass filter to apply to pattern - optional field
    if isfield(settings,'sigmaF')
        sigmaF = settings.sigmaF;
    else
        sigmaF = sind(theta_max)/mean(illum_tbl(:,1)); % Default LPF
    end
    fx = orders_x/px;
    fy = orders_y/py;
    WLPF = exp(-0.5*((fx/sigmaF).^2 + (fy/sigmaF).^2));
    %% 
    if ismember(model_str,{'RCWA','FDTD','FDTD (FBC)'})%strcmp(model_str,'RCWA')
        %% Calculate using Panoramic
        settings_panormic = settings.settings_panormic;
        settings_panormic.px = px;
        settings_panormic.py = py;
        settings_panormic.Dx = settings.Dx;
        settings_panormic.Dy = settings.Dy;
        settings_panormic.orders_x = orders_x;
        settings_panormic.orders_y = orders_y;
        settings_panormic.runType = settings.model_str;
        E = calcPanoramic(model,illum_tbl,settings_panormic); % Note: only compatible with 1D (currently)
    else
        %% Calculate using Fresnel equations and binary pattern
        %% Calcualte pattern and its complement
        if isfield(settings,'P') % Pattern defined in terms of Fourier coefficients
            P = settings.P;
        else  % Binary square-wave pattern
            Dx = settings.Dx;
            Dy = settings.Dy;
            sox = getOrdersRect_background_vec(0,Dx,1,0,orders_x(:));
            sox = sox.';
            soy = getOrdersRect_background_vec(0,Dy,1,0,orders_y(:));

            P = sox.*soy;
        end
        
        if isfield(settings,'Pc')
            Pc = settings.Pc;
        else
            Pc = fft2c(1 - ifft2c(P)*numel(P))/numel(P); % Complement of the pattern
        end
        
        if ~isnan(sigmaF) % Apply low-pass filter to pattern function and its complement
            P = P.*WLPF; % Pattern and its complement
            Pc = Pc.*WLPF;
        end
        %% Calculate Fresnel components
        switch model_str
            case 'TMM a'
                outputs = calcMLAbsEtchFresnel_model(model,illum_tbl);
                r_MLabs = outputs.r(:,1);
                r_MLetch = outputs.r(:,2);
                
                r_MLabs = reshape(r_MLabs,1,1,[]);
                r_MLetch = reshape(r_MLetch,1,1,[]);
                
                E = r_MLetch.*P + r_MLabs.*Pc;
            case 'TMM b'
                rt = fresTFs(model,illum_tbl,orders_x,orders_y,px,py,theta_max);
                
                r_MLabs = rt(:,:,:,1);
                r_MLetch = rt(:,:,:,2);
                
                E = r_MLetch.*P + r_MLabs.*Pc;
%                 r_MLabs = reshape(r_MLabs,1,1,[]);
%                 r_MLetch = reshape(r_MLetch,1,1,[]);
            case 'DblSc'
                rt = fresTFs(model,illum_tbl,orders_x,orders_y,px,py,theta_max);
                
                r_ML = rt(:,:,:,3);
                r_abs = rt(:,:,:,4);
                r_etch = rt(:,:,:,5);
                t_abs = rt(:,:,:,6);
                t_etch = rt(:,:,:,7);
                
                R0 = r_etch.*P + r_abs.*Pc;
                T0 = t_etch.*P + t_abs.*Pc;
%                 T1 = T0 .* r_ML;
                
%                 E = R0 + fft2c( ifft2c( T0*numel(P) ) .* ifft2c( T1 )*numel(P) )/numel(P);
                if ~isfield(settings,'dx')
                    T1 = T0 .* r_ML;
                    E = R0 + fft2c( ifft2c( T0*numel(P) ) .* ifft2c( T1 )*numel(P) )/numel(P);
                else
                    s = size(T0);
                    Nx = round(px/settings.dx/2)*2-1;
                    Nx = max(Nx,max(orders_x));
                    orders_x1 = -Nx:Nx;
                    Ny = round(py/settings.dy/2)*2-1;
                    Ny = max(Ny,max(orders_y));
                    orders_y1 = -Ny:Ny;
                    
                    iix = ismember(orders_x1,orders_x);
                    iiy = ismember(orders_y1,orders_y);
                    
                    s(1) = length(orders_y1);
                    s(2) = length(orders_x1);
                    T0a = zeros(s);
                    r_MLa = zeros(s);
                    T0a(iiy,iix,:,:,:,:) = T0;
                    T0 = T0a;
                    r_MLa(iiy,iix,:,:,:,:) = r_ML;
                    r_ML = r_MLa;
                    T1 = T0 .* r_ML;
                    n = prod(s(1:2));
                    E = fft2c( ifft2c( T0*n ) .* ifft2c( T1*n ) )/n;
                    E = E(iiy,iix,:,:,:,:);
                    E = E + R0;
                end
            case 'DblSc a'
                outputs = calcMLAbsEtchFresnel_model(model,illum_tbl);
                r_ML = fresMLTF(model,illum_tbl,orders_x,orders_y,px,py,theta_max);

%                 r_ML = outputs.r(:,3);
                r_abs = outputs.r(:,4);
                r_etch = outputs.r(:,5);
                t_abs = outputs.t(:,1);
                t_etch = outputs.t(:,2);
                
%                 r_ML = reshape(r_ML,1,1,[]);
                r_abs = reshape(r_abs,1,1,[]);
                r_etch = reshape(r_etch,1,1,[]);
                t_abs = reshape(t_abs,1,1,[]);
                t_etch = reshape(t_etch,1,1,[]);
                
                R0 = r_etch.*P + r_abs.*Pc;
                T0 = t_etch.*P + t_abs.*Pc;
                
                
                if ~isfield(settings,'dx')
                    T1 = T0 .* r_ML;
                    E = R0 + fft2c( ifft2c( T0*numel(P) ) .* ifft2c( T1 )*numel(P) )/numel(P);
                else
                    s = size(T0);
                    Nx = round(px/settings.dx/2)*2-1;
                    Nx = max(Nx,max(orders_x));
                    orders_x1 = -Nx:Nx;
                    Ny = round(py/settings.dy/2)*2-1;
                    Ny = max(Ny,max(orders_y));
                    orders_y1 = -Ny:Ny;
                    
                    iix = ismember(orders_x1,orders_x);
                    iiy = ismember(orders_y1,orders_y);
                    
                    s(1) = length(orders_y1);
                    s(2) = length(orders_x1);
                    T0a = zeros(s);
                    r_MLa = zeros(s);
                    T0a(iiy,iix,:,:,:,:) = T0;
                    T0 = T0a;
                    r_MLa(iiy,iix,:,:,:,:) = r_ML;
                    r_ML = r_MLa;
                    T1 = T0 .* r_ML;
                    n = prod(s(1:2));
                    E = fft2c( ifft2c( T0*n ) .* ifft2c( T1*n ) )/n;
                    E = E(iiy,iix,:,:,:,:);
                    E = E + R0;
                end
        end
    end
end

function E = calcPanoramic(model,illum_tbl,settings_panormic)
    filepath = settings_panormic.filepath;
    blankStr = settings_panormic.blankStr;
    is1D = settings_panormic.is1D;
    Nx_panoramic = settings_panormic.Nx_panoramic;
    dz = settings_panormic.dz; % 0.01
    Nz_plane = settings_panormic.Nz_plane; % 2000
    
    px = settings_panormic.px;
    py = settings_panormic.py;
    
    Dx = settings_panormic.Dx;
    Dy = settings_panormic.Dy;
    
%     runType = settings_panormic.runType;
    
    ind_struct = model.ind_struct;
    x = model.x;
    materials = model.materials;
    lambdaTheta = illum_tbl(:,1:2);
    
    thick = x(ind_struct.thick);
    rough = 0*x(ind_struct.rough);
    composition = x(ind_struct.composition);
    composition = reshape(composition,length(materials),[]);
    [f0f1_elements,rho_nom,amu] = getScatFac(materials,lambdaTheta(:,1));

    n_etch = model.n_etch;
    n_per = model.n_per;
    n_cap = model.n_cap;
    N_ML = model.N_ML;

    out = mkMLAbsEtch(n_cap,n_per,n_etch,thick,amu,rho_nom,rough,f0f1_elements,composition,lambdaTheta);

    layer_thick = [out.absorber.thick;repmat(out.mirror.thick,1,1);thick(end-1:end)];%out.etch.thick(2:3)];
    layer_nk = [out.absorber.nk,repmat(out.mirror.nk,1,1),out.etch.nk(:,2:3)];

    % Label for materials, needs to match EM-Suite template file
    layer_mat = cellfun(@(x) ['M' num2str(x)],num2cell(1:50),'UniformOutput',false);
    layer_mat1 = layer_mat(1:length(layer_thick));
    layer_nk1 = layer_nk;
    layer_nk = [layer_nk1,ones(size(layer_nk1,1),length(layer_mat)-length(layer_mat1))];
    
    %                             layer_mat = cellfun(@(x) ['M' num2str(x)],num2cell(1:length(layer_thick)),'UniformOutput',false);
    if strcmp(settings_panormic.runType,'FDTD (FBC)')
        layer_mat_FBC = cellfun(@(x) ['M' num2str(x) 'FBC'],num2cell(1:4),'UniformOutput',false);
        layer_nk_FBC = out.mirror.nk;
        t_FBC = out.mirror.thick;
        layer_t_FBC = cellfun(@(x) ['t' num2str(x) 'FBC'],num2cell(1:4),'UniformOutput',false);
    else
        layer_mat_FBC = {};
        layer_nk_FBC = [];
        
        layer_t_FBC = {};
        t_FBC = [];
    end
    

%     Nx_panoramic = app.NxEditField.Value;
    
%     dz = 0.01;
%     Nz_plane = 2000;

    simParams = [];

    simParams.Nx = Nx_panoramic;
    simParams.dz = dz;
    simParams.Lz = 0;
    simParams.zIn = simParams.Lz - 2*Nz_plane*dz;
    simParams.zBackward = simParams.Lz - Nz_plane*dz;

    illumParams = [];
    illumParams.lambda = illum_tbl(:,1);
    illumParams.theta = illum_tbl(:,2);
    illumParams.phi = illum_tbl(:,3);


    %             SWA = d(5);
    geometryParams = [];
    geometryParams.t = -sum(layer_thick(1:n_etch)) + layer_thick(end-1);
    if is1D
        %%
        simParams.mag = 1;
        simParams.p = px;

        geometryParams.type = 'box';
        geometryParams.D = Dx;
        geometryParams.z0 = 2*Nz_plane*dz;
        geometryParams.matStr = 'Vacuum';

        geometryParams.box2 = [];
        geometryParams.box2.D = Dx;
        geometryParams.box2.matStr = layer_mat{n_cap+1+n_per+1};
        geometryParams.box2.t = layer_thick(end);
        geometryParams.box2.dz = 0;
    else
        simParams.magx = 1;
        simParams.magy = 1;
        simParams.px = px;
        simParams.py = py;
        
        geometryParams.type = 'box2D';
        
        geometryParams.Dx = Dx;
        geometryParams.Dy = Dy;
        geometryParams.z0 = 2*Nz_plane*dz;
        geometryParams.matStr = 'Vacuum';
        
%         coordsX = geometryParams.coordsX;
%         coordsY = geometryParams.coordsY;
%         Dx = geometryParams.Dx;
%         Dy = geometryParams.Dy;
%         t = geometryParams.t;
%         matStr = geometryParams.matStr;
    end
    %%
    if ~isempty(layer_t_FBC)
        for j = 1:length(layer_t_FBC)
            simParams.(layer_t_FBC{j}) = t_FBC(j);
        end
    end
    %% Connect to API server
    javaaddpath('C:\Program Files\Panoramic\v700\api\MATLAB_6_5.jar');
    connectPanoramic();
    %%
    layer_mat_cap = layer_mat(1:n_cap);
    layer_mat_per = layer_mat(n_cap+(1:n_per));
    %                             layer_mat_cap = layer_mat(1:n_cap+1);
    %                             layer_mat_per = layer_mat(n_cap+1+(1:n_per));

    materialParams = [];
    materialParams.materials = [layer_mat,layer_mat_FBC];
    materialParams.nk = [layer_nk,layer_nk_FBC];
    %                             [materialParams.materials,inds] = unique(layer_mat);
    %                             materialParams.nk = layer_nk(:,inds);
    materialParams.updateNK = false;

    stackParams = [];
    stackParams.Cap = [];
    stackParams.Cap.t = [out.absorber.thick;repmat(out.mirror.thick,N_ML,1)];
    stackParams.Cap.material = [layer_mat_cap,repmat(layer_mat_per,1,N_ML)];

    stackParams.Periodic = [];
    stackParams.Periodic.N = 0;
    stackParams.Periodic.dspace = 1;
    stackParams.Periodic.t = [];
    stackParams.Periodic.material = [];

    % stackParams.Etch = [];
    % stackParams.Etch.t = layer_thick(end-1:end);
    % stackParams.Etch.material = layer_mat_contam;

    geometryParams.stack = stackParams;

    %% Define function for EUV
    struct_EUV = [];

    struct_EUV.simParams = simParams;
    struct_EUV.reflMode = 1;
    struct_EUV.illumParams =  illumParams;
    struct_EUV.materialParams = materialParams;
    struct_EUV.geometryParams = geometryParams;
    struct_EUV.blankStr = blankStr;
%     struct_EUV.blankStr = app.BlankEMSuiteModelFileTextArea.Value{1};
    if isfield(settings_panormic,'runType')
        struct_EUV.runType = settings_panormic.runType;
    else
        struct_EUV.runType = 'RCWA';
    end
    struct_EUV.Nz_plane = Nz_plane;

    nx = simParams.Nx/2;
    orders_x = -nx:nx-1;


    struct_EUV.lambda = illumParams.lambda;
    struct_EUV.theta = illumParams.theta;
    struct_EUV.phi = illumParams.phi;

    if is1D
        struct_EUV.geometryParams.coordsY = [0 10];
        struct_EUV.geometryParams.coordsX = [0 struct_EUV.simParams.mag*struct_EUV.simParams.p];
        struct_EUV.orders = orders_x';
    else
        orders_y = orders_x';
        struct_EUV.geometryParams.coordsX = [0 struct_EUV.simParams.magx*struct_EUV.simParams.px];
        struct_EUV.geometryParams.coordsY = [0 struct_EUV.simParams.magy*struct_EUV.simParams.py];
        struct_EUV.orders_x = orders_x;
        struct_EUV.orders_y = orders_y';
    end

%     filepath = fileparts(app.BlankEMSuiteModelFileTextArea.Value{1});
    struct_EUV.saveLoc = [filepath '\Output.sim'];

    %%
    geometryParams = struct_EUV.geometryParams;
    simParams = struct_EUV.simParams;
    illumParams = struct_EUV.illumParams;
    materialParams = struct_EUV.materialParams;

    Nillum = length(illumParams.lambda);

    if is1D
        eField = zeros(simParams.Nx,Nillum,2,3);
%         scatteredOrders = zeros(length(orders_x),Nillum,2,2);
    else
        eField = zeros(simParams.Nx,simParams.Nx,Nillum,2,3);
%         scatteredOrders = zeros(length(orders_y),length(orders_x),Nillum,2,2);
    end
    %                                 saveLoc = struct_EUV.saveLoc;
    %                                 blankStr = struct_EUV.blankStr;
    %% Load simulation with blank stack
    loadSetup(struct_EUV.blankStr);
    An = getVariableNames(0);
    An = cellfun(@char,cell(An),'UniformOutput',false);
    %                                 loadSetupPanoramic(app,struct_EUV);
    %% Generate the simulation domain
    Nz_plane = struct_EUV.Nz_plane;
    dz = simParams.dz;

    if ~isfield(geometryParams,'stack')
        tML = 0;
    else
        stack = geometryParams.stack;
        tML = sum(stack.Cap.t) + stack.Periodic.N*stack.Periodic.dspace; % Multilayer thickness
    end
    if struct_EUV.reflMode
        geometryParams.z0 = Nz_plane*dz;
        simParams.Lz = round((geometryParams.z0 + max(0,geometryParams.t) + tML + 3*Nz_plane*dz)/dz)*dz; 
        simParams.zIn = simParams.Lz - 2*Nz_plane*dz;
        simParams.zBackward = simParams.Lz - Nz_plane*dz;
    else
        if isfield(geometryParams,'box2')
            delta_z = geometryParams.box2.t + geometryParams.box2.dz;
        else
            delta_z = 0;
        end
        geometryParams.z0 = 2*Nz_plane*dz;
        simParams.Lz = round((geometryParams.z0 + delta_z + max(0,geometryParams.t) + tML + 2*Nz_plane*dz)/dz)*dz; 
        simParams.zIn = simParams.Lz - Nz_plane*dz;
        simParams.zBackward = geometryParams.z0 - Nz_plane*dz;
    end
    if ~isempty(layer_t_FBC)
        simParams.zObs_FBC = geometryParams.z0 - dz;
        simParams.zExc_FBC = geometryParams.z0;
    end
    %% Set simulation parameters
    simFields = fieldnames(simParams);
    for j = 1:length(simFields)
        str = simFields{j};
        if ismember(str,An)
            setVariableValues(str,simParams.(str));
        end
    end

    %% Set illumination parameters
    illumFields = fieldnames(illumParams);
    for j = 1:length(illumFields)
        str = illumFields{j};
        if ismember(str,An)
            setVariableValues(str,illumParams.(str));
        end
    end
    %% Set material parameters
    for j = 1:length(materialParams.materials)
        str = ['n' materialParams.materials{j}];
        if isfield(materialParams,str)
            setVariableValues(str,materialParams.(str))
        else
        if ismember(str,An)
            setVariableValues(str,real(materialParams.nk(:,j)));
        end

        end
        str = ['k' materialParams.materials{j}];
        if isfield(materialParams,str)
            setVariableValues(str,materialParams.(str))
        else
        if ismember(str,An)
            setVariableValues(str,-imag(materialParams.nk(:,j)));
        end

        end
    end

    %% Create the geometry
    %% Set up geometry
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
    elseif isfield(geometryParams,'box2_2D')
        coordsX = geometryParams.box2.coordsX;
        coordsY = geometryParams.box2.coordsY;
        Dx = geometryParams.box2.Dx;
        Dy = geometryParams.box2.Dy;
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
        coordsX_tmp = (coordsX - mean(coordsX))*Dx + mean(coordsX);
        coordsY_tmp = (coordsY - mean(coordsY))*Dy + mean(coordsY);
        coords = [coordsX_tmp coordsY_tmp z0 z0+t];
        coordString = coords2String(coords);
        material = matStr;
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
    end
    %% Save setup
    callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','saveSetup',struct_EUV.saveLoc)
    %% Run
    switch struct_EUV.runType
        case 'RCWA'
            outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingRCWA',1);
        case {'FDTD','FDTD (FBC)'}
            outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingTEMPESTpr2',1);
        case 'TMM'
            outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingStackKirchhoffThick',1);
            
    end
%     %% Collect scattered orders
%     [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
%     [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);
% 
%     %% Organize simulation outputs
%     tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
%     tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
%     cgrp = unique(tableTE(:,strcmp(header,'cgrp')));
%     %% Pull out the orders 
%     for iCGRP = 1:length(cgrp)
%         col_ind = 3;
%         rows = tableTE(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
%         if length(header) > 3
%             col_ind = strcmp(header,'Data');
%         end
%         %%
%         orders_tmp = tableTE(rows,strcmp(header,'m'));
%         %                                     orders_tmp = reshape(orders_tmp,[],Nillum);
%         orders_tmp = orders_tmp - 0.5;
%         ETE = tableTE(rows,col_ind);
%         ETM = tableTM(rows,col_ind);
%         ind0 = 0;
%         for j = 1:Nillum
%             ntmp = find(abs(diff(orders_tmp(ind0+1:end))) > 2,1);
%             inds = ind0+(1:ntmp);
%             [orders_u,iA,iC] = intersect(orders_tmp(inds),orders_x);
%             %% TE Output
%             scatteredOrders(iC,j,iCGRP,1) = ETE(inds(iA));
%             %% TM Output
%             scatteredOrders(iC,j,iCGRP,2) = ETM(inds(iA));
%         end
%     end

    %% Collect E-fields
    [tableEx,~] = getEMOutputTable(outh,'Ex',true);
    [tableEy,~] = getEMOutputTable(outh,'Ey',true);
    [tableEz,header] = getEMOutputTable(outh,'Ez',true);
    cgrp = unique(tableEz(:,strcmp(header,'cgrp')));

    if is1D
        xx = unique(tableEx(:,strcmp(header,'x')));
        for iCGRP = 1:length(cgrp)
            col_ind = 3;
            rows = tableEx(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            if length(header) > 3
                col_ind = strcmp(header,'Data');
            end
            %%


            Ex = reshape(tableEx(rows,col_ind),[],Nillum);
            Ey = reshape(tableEy(rows,col_ind),[],Nillum);
            Ez = reshape(tableEz(rows,col_ind),[],Nillum);

            eField(:,:,iCGRP,1) = Ex;
            eField(:,:,iCGRP,2) = Ey;
            eField(:,:,iCGRP,3) = Ez;
        end
    else
        xx = unique(tableEx(:,strcmp(header,'x')));
        yy = unique(tableEx(:,strcmp(header,'y')));
        for iCGRP = 1:length(cgrp)
            col_ind = 3;
            rows = tableEx(:,strcmp(header,'cgrp')) == cgrp(iCGRP);
            if length(header) > 3
                col_ind = strcmp(header,'Data');
            end
            %%


            Ex = reshape(tableEx(rows,col_ind),length(yy),length(xx),Nillum);
            Ey = reshape(tableEy(rows,col_ind),length(yy),length(xx),Nillum);
            Ez = reshape(tableEz(rows,col_ind),length(yy),length(xx),Nillum);

            eField(:,:,:,iCGRP,1) = Ex;
            eField(:,:,:,iCGRP,2) = Ey;
            eField(:,:,:,iCGRP,3) = Ez;
        end
    end
    %%
    if is1D
        N_full = length(xx)/2;
        orders_full = -N_full:N_full-1;

        lambda = illum_tbl(:,1);
        theta = illum_tbl(:,2);
        phi = illum_tbl(:,3);
        phaseRamp = 2*pi*xx./lambda.'.*sind(theta.').*cosd(phi.');
        %                                 phaseRamp = 2*pi*xx/lambda(i)*sind(theta(i))*cosd(phi(i));
        phaseRampX = phaseRamp;
        phaseRampY = phaseRamp;

        Ex = eField(:,:,1,1).*exp(1i*phaseRampX);
        Ey = -eField(:,:,1,2).*exp(1i*phaseRampY);
        E =  Ex.*sind(phi.') + Ey.*cosd(phi.');
        E_FT = fftc(E,1)/length(E);
        ii = ismember(orders_full,orders_x);

        soRCWA = E_FT(ii,:);

        %% Correct phase for propagation of illumination and scattered wave
        p = simParams.p*simParams.mag;
        k0 = 2*pi./lambda.';
        ky = k0.*sind(theta.').*sind(phi.');
        kx = k0.*sind(theta.').*cosd(phi.') - 2*pi*orders_x(:)./p;

        %                                 ii = abs(kx.^2 + ky.^2) <= k0.^2;
        kz = real(sqrt(k0.^2 - kx.^2 - ky.^2));
        %                                 [kx,ky,kz] = getKVec_1D(orders_x(:),geometryParams.coordsX(2),illumParams.lambda.',illumParams.theta.',illumParams.phi.');
        if struct_EUV.reflMode
            zT = geometryParams.z0 + max(0,geometryParams.t) + tML;
            z1 = simParams.zIn - zT;
            z2 = simParams.zBackward - zT;
            soRCWA = soRCWA.*exp(1i*kz(orders_x(:) == 0,:)*z1).*exp(1i*kz*z2);
        else
            zT = geometryParams.z0 + max(0,geometryParams.t) + tML;
            zB = geometryParams.z0 + tML;
            z1 = simParams.zIn - zT;
            z2 = zB - simParams.zBackward;
            soRCWA = soRCWA.*exp(1i*kz(orders_x(:) == 0,:)*z1).*exp(1i*kz*z2);
        end
    else
        Nx_full = length(xx)/2;
        orders_x_full = -Nx_full:Nx_full-1;
        Ny_full = length(yy)/2;
        orders_y_full = -Ny_full:Ny_full-1;

        lambda = illum_tbl(:,1);
        lambda = reshape(lambda,1,1,[]);
        theta = illum_tbl(:,2);
        theta = reshape(theta,1,1,[]);
        phi = illum_tbl(:,3);
        phi = reshape(phi,1,1,[]);
        phaseRamp = 2*pi./lambda.*sind(theta).*(xx.'.*cosd(phi) + yy.*sind(phi));
        %                                 phaseRamp = 2*pi*xx/lambda(i)*sind(theta(i))*cosd(phi(i));
        phaseRampX = phaseRamp;
        phaseRampY = phaseRamp;

        Ex = eField(:,:,1,1).*exp(1i*phaseRampX);
        Ey = -eField(:,:,1,2).*exp(1i*phaseRampY);
        E =  Ex.*sind(phi) + Ey.*cosd(phi);
        E_FT = fft2c(E)/numel(E(:,:,1));
        ii_x = ismember(orders_x_full,orders_x);
        ii_y = ismember(orders_y_full,orders_y);

        soRCWA = E_FT(ii_y,ii_x,:);

        %% Correct phase for propagation of illumination and scattered wave
        px = simParams.px*simParams.magx;
        py = simParams.py*simParams.magy;
        k0 = 2*pi./lambda;
%         ky = k0.*sind(theta).*sind(phi);
        kx = k0.*sind(theta).*cosd(phi) - 2*pi*orders_x./px;
        ky = k0.*sind(theta).*sind(phi) - 2*pi*orders_y./py;

        %%                                 ii = abs(kx.^2 + ky.^2) <= k0.^2;
        kz = real(sqrt(k0.^2 - kx.^2 - ky.^2));
        %                                 [kx,ky,kz] = getKVec_1D(orders_x(:),geometryParams.coordsX(2),illumParams.lambda.',illumParams.theta.',illumParams.phi.');
        if struct_EUV.reflMode
            zT = geometryParams.z0 + max(0,geometryParams.t) + tML;
            z1 = simParams.zIn - zT;
            z2 = simParams.zBackward - zT;
            soRCWA = soRCWA.*exp(1i*kz(orders_y == 0,orders_x==0,:)*z1).*exp(1i*kz*z2);
        else
            zT = geometryParams.z0 + max(0,geometryParams.t) + tML;
            zB = geometryParams.z0 + tML;
            z1 = simParams.zIn - zT;
            z2 = zB - simParams.zBackward;
            soRCWA = soRCWA.*exp(1i*kz(orders_y == 0,orders_x,:)*z1).*exp(1i*kz*z2);
        end
    end
    %% Delete results in Panormaic for memory management
    destroyAll();
    %%
    if is1D
        E = circshift(flipud(squeeze(soRCWA)),mod(size(soRCWA,1)+1,2));
        if length(orders_full)>length(settings_panormic.orders_x)
            E = E(ismember(orders_full,settings_panormic.orders_x),:);
        else
            Etmp = E;
            s = size(E);
            s(1) = length(settings_panormic.orders_y);
            s(2) = length(settings_panormic.orders_x);
            E = zeros(s);
            E(:,ismember(settings_panormic.orders_x,orders_full),:,:,:) = Etmp;
        end
        E = reshape(E,length(settings_panormic.orders_y),length(settings_panormic.orders_x),[]);
    else
        E = soRCWA;
%         E = rot90(soRCWA,1);
        if length(orders_x_full)>length(settings_panormic.orders_x)
%         E = circshift(flipud(fliplr(soRCWA)),mod(size(soRCWA(:,:,1))+1,2));
            E = E(ismember(orders_y_full,settings_panormic.orders_y),ismember(orders_x_full,settings_panormic.orders_x),:);
        else
            Etmp = E;
            s = size(E);
            s(1) = length(settings_panormic.orders_y);
            s(2) = length(settings_panormic.orders_x);
            E = zeros(s);
            E(ismember(settings_panormic.orders_y,orders_y_full),ismember(settings_panormic.orders_x,orders_x_full),:,:,:) = Etmp;
        end
        E = rot90(E,1);
        E = reshape(E,length(settings_panormic.orders_y),length(settings_panormic.orders_x),[]);
    end
end