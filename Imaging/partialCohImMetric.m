function [metric,outputs] = partialCohImMetric(scatter_settings,im_settings,exposure_settings)
    outputs = [];
    %% Imaging settings
    mag_xy = im_settings.mag_xy;
%     NA_mask = im_settings.NA_mask;
    NA_im = im_settings.NA_im;
    NA_mask = im_settings.NA_im./mag_xy;
    zz = im_settings.zz;
    theta0 = im_settings.theta0;
    uxy0 = [0,sind(theta0)]./NA_mask;
    
    source = im_settings.source;
%     nSV = im_settings.nSV;
    
    %% Near-field scattering settings
    settings = scatter_settings.settings;
    illum_tbl = scatter_settings.illum_tbl;
    settings.illum_tbl = illum_tbl;
    orders_x = settings.orders_x;
    orders_y = settings.orders_y;
    px = settings.px*mag_xy(1);
    py = settings.py*mag_xy(2);
    %% Exposure settings
    if isfield(exposure_settings,'area_target')
        area_target = exposure_settings.area_target;
    end
    exposure_lattitude = exposure_settings.exposure_lattitude;
    area_lattitude = exposure_settings.area_lattitude;
    if isfield(exposure_settings,'metric')
        metric = exposure_settings.metric;
    else
        metric = 'Area';
    end
    if strcmp(metric,'Target correlation')
        P_target = exposure_settings.P_target;
    elseif strcmp(metric,'Print no print')
        P_print = exposure_settings.P_print;
        P_dontprint = exposure_settings.P_dontprint;
    end
    if isfield(exposure_settings,'st')
        st = exposure_settings.st;
    else
        st = strel('disk',2);
    end
    if isfield(exposure_settings,'metric1')
        metric1 = exposure_settings.metric1;
        if isfield(exposure_settings,'alpha1')
            alpha1 = exposure_settings.alpha1;
        end
    else
        metric1 = '';
    end
    if isfield(exposure_settings,'Wpw')
        Wpw = exposure_settings.Wpw;
    else
        Wpw = ones(length(zz),length(exposure_lattitude));
    end
    Wpw = reshape(Wpw,1,1,length(zz),length(exposure_lattitude));
    Wpw = Wpw/mean(Wpw(:));
    %% Far-field focus stack
    if isfield(im_settings,'z0')
        z0 = im_settings.z0;
    else
        z0 = 0;
    end
    if isfield(im_settings,'z0_fit')
        z0_fit = im_settings.z0_fit;
    else
        z0_fit = false;
    end
    if isfield(im_settings,'LBz')
        LBz = im_settings.LBz;
    else
        LBz = -100;
    end
    if isfield(im_settings,'UBz')
        UBz = im_settings.UBz;
    else
        UBz = 100;
    end
    if isfield(im_settings,'I')
        I = im_settings.I;
        outputs.I_PC = I;
    elseif isfield(im_settings,'I_PC')
        I_PC = im_settings.I_PC;
        
        if false
            if isfield(im_settings,'w')
                w = im_settings.w;
            else
                if isfield(im_settings,'V')
                    V = im_settings.V;
                else
                    V = eye(numel(source));
                end
                w = V'*source(:);%/length(source);
            end
            w = w/norm(w);
    %         I_PC = zeros(size(E,1),size(E,2),length(zz),nSV);
    %         for j = 1:size(U,2)
    %             tmp = U(:,j)*S(j,j);
    %             tmp = reshape(tmp,size(I_PC(:,:,:,j)));
    %             I_PC(:,:,:,j) = abs(ifft2c(tmp*length(orders_y)*length(orders_x))).^2;
    %         end
            I = sum(reshape(abs(w),1,1,1,[]).*I_PC,4);
            I = reshape(I,length(orders_y),length(orders_x),[]);
        else
            I = mean(source.*I_PC,3)./mean(source,3);
            I = reshape(I,length(orders_y),length(orders_x),[]);
        end
        
        if isfield(exposure_settings,'I0')
            I0 = exposure_settings.I0;
        else
            I0 = mean(I,'all');
        end
        
        
        
    else
        %% Near-field scattering for each source point
        if isfield(scatter_settings,'E')
            E = scatter_settings.E;
        else
            E = nearFieldDoubleScattering(settings);
        end
        outputs.E = E;
%         %% Illumination grid in mask space
%         illum_tbl2 = getDiffractionIllum2D(illum_tbl,orders_x,orders_y,px,py);
%         kxyz2 = illum2kxyz(illum_tbl2);
%         uxy2 = kxyz2(:,1:2)./(2*pi./illum_tbl2(:,1))./NA_mask - uxy0;
%         kp = rssq(uxy2,2) <= 1;

%         H = reshape(kp,size(illum_tbl,1),length(orders_y),length(orders_x));
%         H = permute(H,[2 3 1]);
        %% Convert to image space
        px = settings.px;
        py = settings.py;

        kxyz = illum2kxyz(illum_tbl);
        k0 = 2*pi./illum_tbl(:,1);
        uxy = kxyz(:,1:2)./k0./NA_mask - uxy0;
        
        illum_tbl_im = [illum_tbl(:,1),real(asind(rssq(uxy*NA_im,2))),atan2d(uxy(:,2),uxy(:,1))];
        
        illum_tbl2 = getDiffractionIllum2D(illum_tbl_im,orders_x,orders_y,px,py);
        
        kxyz2 = illum2kxyz(illum_tbl2);
        uxy2 = kxyz2(:,1:2)./(2*pi./illum_tbl2(:,1));
%         uxy2 = uxy.*NA_im;
        
        illum_tbl_im2 = illum_tbl2;

        kp = rssq(uxy2,2) <= NA_im;

        H = reshape(kp,size(illum_tbl,1),length(orders_y),length(orders_x));
        H = permute(H,[2 3 1]);
        
%         kxyz(:,1:2) = k0.*uxy2;
%         kxyz(:,3) = real( sqrt( k0.^2 - sum(kxyz(:,1:2).^2,2)));
%         illum_tbl_im = kxyz2illum(kxyz);
%         illum_tbl_im2 = getDiffractionIllum2D(illum_tbl_im,orders_x,orders_y,px,py);
    %     %%
    %     
    %     kxyz = kxyz - mean(kxyz);
    % %     kxyz = kxyz2;
    %     kxyz(:,1:2) = 2*pi*NA_im*uxy2;
    %     kxyz(:,3) = real( sqrt( (2*pi./illum_tbl2(:,1)).^2 - sum(kxyz(:,1:2).^2,2)));
    %     illum_tbl_im = kxyz2illum(kxyz);
    %     illum_tbl_im2 = getDiffractionIllum2D(illum_tbl_im,orders_x,orders_y,px,py);
        kxyz_im2 = illum2kxyz(illum_tbl_im2);
        %% Compute defocused and bandlimitted far-field for each source point
        fHz = @(z0) H.*exp(1i*reshape(kxyz_im2(:,3),length(orders_y),length(orders_x),[]).*reshape(zz+z0,1,1,1,[]));
        fImPC = @(E,Hz) abs(ifft2c(E.*Hz*length(orders_y)*length(orders_x))).^2;
        fI = @(I_PC,source) reshape(mean(source.*I_PC,3)./mean(source,3),length(orders_y),length(orders_x),[]);
        fz02I = @(z0) fI(fImPC(E,fHz(z0)),source);
        
        if strcmp(metric,'Target correlation')
            fFEMerr = @(thr,z0) sqrt(mean(((fz02I(z0)>=thr.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
    %         fFEMerr = @(x) abs(mean((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
        elseif strcmp(metric,'Area')
            fFEMerr = @(thr,z0) abs(mean((fz02I(z0)>=thr.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
        elseif strcmp(metric,'Print no print')
            fFEMerr = @(thr,z0) err_printnoprint(fz02I(z0),thr,exposure_lattitude,P_print,P_dontprint,orders_x,orders_y,Wpw);
    %         fFEMerr = @(x) mean(abs((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
        end

    %     st = strel('disk',2);
        fPW = @(x) imdilate( imerode( fFEMerr(x(1),x(2))<= area_lattitude,st),st);
        fPWerr = @(x) -mean( fPW(x) ,'all');

        if strcmp(metric1,'EPE')
            ferr = @(x) fFEMerr(x(1),x(2)) + alpha1*err_EPE(fz02I(z0),x,exposure_lattitude,orders_x/length(orders_x),orders_y/length(orders_y));
        else
            ferr = @(x) fFEMerr(x(1),x(2));
        end
        
        if ~z0_fit
            I_PC = fImPC(E,fHz(z0));
            I = fz02I(0);
%             Hz = hHz(0);
%             I_PC = fImPC(E,Hz);
%             I_PC = abs(ifft2c(E.*Hz*length(orders_y)*length(orders_x))).^2;
%             I = mean(source.*I_PC,3)./mean(source,3);
%             I = reshape(I,length(orders_y),length(orders_x),[]);
            if strcmp(metric,'Target correlation')
                fFEMerr = @(x) sqrt(mean(((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
            %         fFEMerr = @(x) abs(mean((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
            elseif strcmp(metric,'Area')
                fFEMerr = @(x) abs(mean((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
            elseif strcmp(metric,'Print no print')
                fFEMerr = @(x) err_printnoprint(I,x,exposure_lattitude,P_print,P_dontprint,orders_x,orders_y,Wpw);
            %         fFEMerr = @(x) mean(abs((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
            end

            %     st = strel('disk',2);
            fPW = @(x) imdilate( imerode( fFEMerr(x)<= area_lattitude,st),st);
            fPWerr = @(x) -mean( fPW(x) ,'all');

            if strcmp(metric1,'EPE')
                ferr = @(x) fFEMerr(x) + alpha1*err_EPE(I,x,exposure_lattitude,orders_x/length(orders_x),orders_y/length(orders_y));
            else
                ferr = @(x) fFEMerr(x);
            end
            %
            %     fPW = @(x) mean(fFEMerr(x)<=area_lattitude,'all');
            II = linspace(min(I(:))+(max(I(:))-min(I(:)))*0.1,max(I(:))-(max(I(:))-min(I(:)))*0.1,101);
            cc = 0*II;
            for i = 1:length(II)
                cc(i) = mean(ferr(II(i)),'all');
            %         cc(i) = fPW(II(i));
            end

            [~,i] = min(cc);
            I0 = II(i);
            thr_best = I0;
            %     W = cc <= min(cc)/0.9;
            %     plot(cc)
            %     W = cc >= max(cc)*0.9;
            %     I0 = mean(W.*II)/mean(W);
            FEM = squeeze(ferr(I0));
            W = double(FEM <= area_lattitude);
            %%
            if sum(W) > 0
                I0 = I0*mean(W.*(1+exposure_lattitude),[1 2])/mean(W,[1 2]);
            end
        else
            I = fz02I(0);
            I0 = mean(I,'all');
            x0 = [I0,z0];
            LB = [min(I,[],'all'),LBz];
            UB = [max(I,[],'all'),UBz];
            ncycles = 4;
            niter = 5;
            fac = 0.75;
            alpha = 0.75;
            step0 = 0.1;
            
            x_fit = gss_coord_descentEWMA(x0,@(x) mean(Wpw.*ferr(x)),LB,UB,ncycles,niter,fac,alpha,step0);
            
%             I_PC = fz02I(x_fit(1));
            
            I0 = x_fit(1);
            z0 = x_fit(2);
%             FEM = squeeze(ferr(x_fit));
            I_PC = fImPC(E,fHz(z0));
        end
        
    %         E_all_FF_z(:,:,:,:,iP) = E(:,:,:,iP).*Hz;
%         %% Compute weighted sum of source points
%         I_PC = abs(ifft2c(E_FF_z*length(orders_y)*length(orders_x))).^2;
%         I = mean(source.*I_PC,3)./mean(source,3);
%         I = reshape(I,length(orders_y),length(orders_x),[]);
%         outputs.I_PC = I_PC;
%         I = fImPC(E,fHz(z0));
        I = fz02I(z0);
%         PW = fPW(x_fit);
    end
    outputs.I = I;
    outputs.I_PC = I_PC;
    %%
    if strcmp(metric,'Target correlation')
        fFEMerr = @(x) sqrt(mean(((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
    %         fFEMerr = @(x) abs(mean((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
    elseif strcmp(metric,'Area')
        fFEMerr = @(x) abs(mean((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))),[1 2])-area_target);
    elseif strcmp(metric,'Print no print')
        fFEMerr = @(x) err_printnoprint(I,x,exposure_lattitude,P_print,P_dontprint,orders_x,orders_y,Wpw);
    %         fFEMerr = @(x) mean(abs((I>=x.*(1+reshape(exposure_lattitude,1,1,1,[]))) - P_target).^2,[1 2]));
    end

    %     st = strel('disk',2);
    fPW = @(x) imdilate( imerode( fFEMerr(x)<= area_lattitude,st),st);
    fPWerr = @(x) -mean( fPW(x) ,'all');

    if strcmp(metric1,'EPE')
        ferr = @(x) fFEMerr(x) + alpha1*err_EPE(I,x,exposure_lattitude,orders_x/length(orders_x),orders_y/length(orders_y));
    else
        ferr = @(x) fFEMerr(x);
    end
    %% Determine best exposure dose
%     I0;
    thr_best = gss_coord_descentEWMA(I0,@(x) mean(Wpw.*ferr(x),'all'),min(I(:)),max(I(:)),10,10,0.9,0.9,0.1);
%     thr_best = I0;
%     thr_best = gss(min(I(:))+(max(I(:))-min(I(:)))*0.1,max(I(:))-(max(I(:))-min(I(:)))*0.1,@(x) -fPW(x),10);
%     thr_best = gss_coord_descentEWMA(I0,@(x) -fPW(x),min(I(:)),max(I(:)),10,10,0.9,0.9,0.1);
    outputs.thr_best = thr_best;
    FEM = squeeze(ferr(thr_best));
    %% Calculate focus and exposure lattitude
    area_error = mean(FEM,'all');
    outputs.area_error = area_error;
    outputs.FEM = FEM;
    outputs.z0 = z0;
    PW = fPW(thr_best);
    outputs.PW = PW;
    if isfield(exposure_settings,'metric_PWerr')
        metric_PWerr = exposure_settings.metric_PWerr;
    else
        metric_PWerr = 'Mean PW';
    end
    
    if strcmp(metric_PWerr,'Mean PW')
        metric = mean(PW.*FEM,'all')/mean(PW,'all') - mean(PW,'all'); % Metric: Mean over range of process window
        if isnan(metric)
            metric = min(FEM(:)); % Metric: minimum value achieved (if PW was empty)
        end 
    elseif strcmp(metric_PWerr,'Mean err')
        metric = mean(squeeze(Wpw).*FEM,'all');
    end
%     metric = mean(PW.*FEM,'all')/mean(PW,'all') - mean(PW,'all'); % Metric: Mean over range of process window
       
    if isfield(exposure_settings,'I0metric')
        metric = exposure_settings.I0metric(metric,thr_best);
    end
%     metric = fPWerr(thr_best);
%         print_ok = area_error <= area_lattitude;
%     end
end

function err = err_printnoprint(I,thr,exposure_lattitude,P_print,P_dontprint,orders_x,orders_y,Wpw)
    %%
    Ibin = double(I>=thr.*(1+reshape(exposure_lattitude,1,1,1,[])));
    I = mean(Ibin,[3 4]);
    N = numel(Ibin(:,:,1,1));
    tmp = real(ifft2c(fft2c(I/N).*conj(fft2c(P_print/N))))*N;
    [~,i0] = max(tmp(:));
    [i0y,i0x] = ind2sub(size(tmp),i0);
    Ibin = circshift(Ibin,-[orders_y(i0y),orders_x(i0x)]);
    err = mean(abs(Ibin - P_print).*P_print,[1 2])/mean(P_print(:));
    err = err + mean(abs(1 - Ibin - P_dontprint).*P_dontprint,[1 2])/mean(P_dontprint(:));
    err = reshape(err,size(err,3),size(err,4));
end

function err = err_EPE(I,thr,exposure_lattitude,x,y)
    %%
    Ibin = double(I>=thr.*(1+reshape(exposure_lattitude,1,1,1,[])));
    if mean(Ibin,'all') > 0
    Ibin = Ibin/mean(Ibin,'all');
    end
    X = mean(Ibin.*x,[1 2]);
    Y = mean(Ibin.*y,[1 2]);
    
    COMx = mean(X,'all');
    COMy = mean(Y,'all');
    
    err = sqrt(mean((X - COMx).^2 + (Y - COMy).^2,'all'));
end


