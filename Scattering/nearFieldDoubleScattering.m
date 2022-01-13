function E = nearFieldDoubleScattering(settings)
    %% Define the thin-film model (absorber, multilayer, etch) - required field
    model = settings.model; % Thin-film model
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
    %% Bandlimiting aperture (optional)
    if isfield(settings,'bandlim')
        theta_xy0 = settings.bandlim.theta_xy0;
        NA = settings.bandlim.NA;
    else
        theta_xy0 = [0 0];
        NA = 1;
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
    rt = fresTFs(model,illum_tbl,orders_x,orders_y,px,py,theta_max,theta_xy0,NA);

    r_ML = rt(:,:,:,3);
    r_abs = rt(:,:,:,4);
    r_etch = rt(:,:,:,5);
    t_abs = rt(:,:,:,6);
    t_etch = rt(:,:,:,7);

    %% Calculate pattern components
    R0 = r_etch.*P + r_abs.*Pc;
    T0 = t_etch.*P + t_abs.*Pc;
    T1 = T0 .* r_ML;

    %% Calculate final field
    E = R0 + fft2c( ifft2c( T0*numel(P) ) .* ifft2c( T1 )*numel(P) )/numel(P);
end