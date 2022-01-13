function [E,X] = solvePhaseLift1D(struct)
    ims = struct.I;
    ims = ims - min(ims);
    ims = ims ./ mean(ims);
    
    y = struct.y1;
    p = struct.p;
    orders = struct.orders1;
    phi_DC = struct.phi_DC;
    zz = struct.zz;
    lambda_scan = struct.lambda_scan;
    %% Read settings from UI
    d = app.UITable3.Data;

    D_nom = d(1);
    r_tbl = d(2)*exp(1i*rad2deg(d(3)));

    N_iter = d(4);
    step_mag = d(5);
    alpha_lr = d(6); % X low-rank prior
    step_decay = d(7);
    N_GD = d(8);
    spinv = d(9);
    sigma_z = d(10);

    Wz = exp(-0.5*(zz(:)'/sigma_z).^2);
    W = repmat(Wz,length(orders),length(phi_DC));
    Weight = vec(W);
    %             Weight = reshape( W,size(W,1),[]);
    %             ims = mean_pad(ims,W);

    NA = app.NAEditField.Value;
    theta0 = app.thetaEditField.Value;

    %% Initialization
    switch app.InitmodeListBox.Value
        case 'Sq wave, specified r'
            r_refl = r_tbl;
        case 'Sq wave, reflectometry'
            s = app.OutputdataTextArea_5.Value{1};
            if ~exist(s,'file')
                return
            end
            tmp = load(s);
            if isfield(tmp,'model')
                model = tmp.model;
            else
                return
            end

            [r1,r2,r_ML,r_abs,t_abs,r_etch,t_unetch] = model.fFres_pol(model.x,'s');

            r = r1./r2;
            S_refl = scatteredInterpolant(model.lambdaTheta(:,1),model.lambdaTheta(:,2),r,'linear');
    end

    %% Coordinate grid
    %             lambda = reshape(lambda_scan,1,1,[]);
    Emat = zeros(length(orders),length(lambda_scan));
    for iL = 1:length(lambda_scan)
        lambda = lambda_scan(iL);

        fc = NA./lambda;
        fx = orders/p;
        ux = fx./fc;
    %                 ux = repmat(ux,1,length(zz),1);
        in_pup = abs(ux) <= 1;
        k0 = 2*pi./lambda;
        kx = k0.*real(sind(asind(sind(theta0)+orders.*lambda/p)));
        kx = kx - kx(orders==0,:,:);
        kx = kx .*in_pup;
    %             kx(~in_pup) = 0;
        kz = k0.*real(cosd(asind(sind(theta0)+orders.*lambda/p)));
        kz = kz - kz(orders==0,:,:);
        kz = kz - orders./2.*(kz(orders==1,:,:)-kz(orders==-1,:,:));
        kz = kz .*in_pup;
    %             kz(~in_pup) = 0;
        H = in_pup.*exp(1i.*kz.*reshape(zz,1,[]));
        H = repmat(H,1,1,length(phi_DC));
        H(orders==0,:,:) = H(orders==0,:,:).*exp(1i*reshape(phi_DC,1,1,[]));
        H = H.*exp(-1i*angle(H(orders==0,:,:)));
        %% PhaseLift matrix
        [B,im_modes,H_modes,F_modes,Basis_mat] = formPhaseLift1D(orders,reshape(H,size(H,1),[]),orders);
        %% Invert the matrix (regularized pseudo-inverse)
        [U,S,V] = svd(B,'econ');
        s = diag(S);
        s(s <= s(1)*spinv) = 0;
        Bpinv = V*diag(s)*U';
        %% Initialization
        switch app.InitmodeListBox.Value
            case 'Sq wave, specified r'
                r_refl = r_tbl;
            case 'Sq wave, reflectometry'
                r_refl = S_refl(lambda,theta0);
        end
        so = getOrdersRect_background_vec(0.5,1-D_nom,r_refl,1,orders(im_modes)');
        so = so*exp(-1i*angle(so(orders(im_modes)==0)));
        so = so/norm(so);
        E0 = so;
        X0 = so*so';
        X0 = X0 / trace(X0);
        %% Data
        b = vec(ims(:,:,iL,:));
        b = mean_pad(b,Weight);
        b = b/mean(b)*trace(X0);
        %%
        X = X0;
        cc = zeros(N_iter,N_GD);
        for i = 1:N_iter
            for j = 1:N_GD
                y = B*vec(X);
                y = mean_pad(y,Weight);
                res = y - b;
                res = res.*Weight;

                cc(i,j) = rms(res)/rms(b);

                semilogy(app.UIAxes2_7,vec(cc'))
                title(app.UIAxes2_7,'MSE')
                drawnow()
                step = Bpinv*res;
                step = step/norm(step)*step_mag*step_decay^(i-1);
                step = reshape(step,size(X));
                X = X - step;
            end
            %% Low-rank
            [U,S,V] = svd(X,'econ');
            s = diag(S);
            s = softThresh(s,alpha_lr*s(1));
            S = diag(s);
            X = U*S*U';
            %% Extract rank 1 component (recovered field)
            E = U(:,1);
            E = E/exp(1i*angle(E(orders(im_modes)==0)));
            E = E/norm(E);
            %% Trace
            X = X/trace(X);
            %% Update plots
            imagesc(app.UIAxes2_5,reshape(b,length(orders),[]))
            app.UIAxes2_5.XLim = [-inf inf];
            app.UIAxes2_5.YLim = [-inf inf];
            app.UIAxes2_5;
            title(app.UIAxes2_5,'Raw data')
            bhat = real(B*vec(E*E'));
            bhat = mean_pad(bhat,Weight);

            imagesc(app.UIAxes2_6,reshape(bhat,length(orders),[]))
            app.UIAxes2_6.XLim  = [-inf inf];
            app.UIAxes2_6.YLim = [-inf inf];
            title(app.UIAxes2_6,'Reconstruction')
            app.UIAxes2_6;

            plot(app.UIAxes2_8,E,'o-','LineWidth',3)
            title(app.UIAxes2_8,'E')
            app.UIAxes2_8;

            drawnow()
        end
        Emat(in_pup,iL) = E;

        if exist(app.OutputdataTextArea.Value{1},'file')
            save(app.OutputdataTextArea.Value{1},'Emat','-append')
        else
            save(app.OutputdataTextArea.Value{1},'Emat')
        end
        updatedphiplot(app)
    end
end