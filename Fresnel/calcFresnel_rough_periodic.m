function [r,t,M_out] = calcFresnel_rough_periodic(n_ML,n_cap,layer_thick,layer_nk,layer_rough,lambdaTheta,pol,n_per)
    if nargin < 8
        n_per = size(layer_thick,1) - n_cap;
    end
%     if ~isrow(layer_thick)
%         layer_thick = layer_thick.';
%     end
    if size(layer_thick,1) == size(layer_nk,2)
        one = ones(size(layer_nk,1),1,size(layer_nk,3));
        layer_nk = [one layer_nk one];
    end
    
    %% Pre-calculate factors 
    nkMat2 = layer_nk.^2;
    alpha2 = nkMat2(:,1).*sind(lambdaTheta(:,2)).^2;
    eta0 = sqrt(nkMat2 - alpha2);
    if pol == 's'
        eta = eta0;
    elseif pol == 'p'
        eta = nkMat2./eta0;
    end
    k0 = 2*pi * 1./lambdaTheta(:,1);
    kz = -k0 .* eta0;
%     kz = 2*pi * 1./lambdaTheta(:,1) .* eta0;
    phi = reshape(layer_thick,size(kz(1,2:end-1,:))) .* kz(:,2:end-1,:);
    
    %%
    kz1 = -k0 .* eta; % effective kz in medium (polarization-dependent)
    %% Calculate transfer matrix
    p =  (kz1(:,2:end,:) + kz1(:,1:end-1,:))./(2*kz1(:,1:end-1,:)); % Symmetric interface term
    m = -(kz1(:,2:end,:) - kz1(:,1:end-1,:))./(2*kz1(:,1:end-1,:)); % Asymmetric interface term

    z_rms = reshape(layer_rough,[1 size(layer_rough)]);
    R = zeros(2,2,size(layer_thick,1)+1,size(lambdaTheta,1),size(layer_thick,2));
    T = zeros(2,2,size(layer_thick,1),size(lambdaTheta,1),size(layer_thick,2));

    %% Refraction matrix (with roughness)
    diff_kz = kz(:,2:end,:) - kz(:,1:end-1,:);
    sum_kz = kz(:,2:end,:) + kz(:,1:end-1,:);
    a = [1 1 size(permute(p, [2 1 3] ))];
    p = p .* exp(-((diff_kz).^2).*(z_rms.^2)/2);
    p = reshape( permute( (p ), [2 1 3] ) , a );
    m = m .* exp(-((sum_kz).^2).*(z_rms.^2)/2);
    m = reshape( permute( (m ), [2 1 3] ) , a );
    R(1,1,:,:,:) = p;
    R(1,2,:,:,:) = m;
    R(2,1,:,:,:) = m;
    R(2,2,:,:,:) = p;

    %% Translation matrix (i.e. bulk propagation)   
    a = [1 1 size(permute(phi, [2 1 3] ))];
    T(1,1,:,:,:) = reshape(  permute( exp(-1i*phi), [2 1 3] ) , a );
    T(2,2,:,:,:) = reshape(  permute( exp(1i*phi), [2 1 3] ) , a );

    %% Calculate transfer matrix (sequential matrix multiplication)
    M_out = reshape(R(:,:,1,:,:),2,2,[]);
    % Capping layers
    for k = 1:n_cap
        M_out = matMul2x2( M_out , T(:,:,k,:,:) );
        M_out = matMul2x2( M_out , R(:,:,k + 1,:,:) );
    end
    % Periodic layers
    if n_ML > 0 
        M = repmat(eye(2,2),size(R(1,1,1,:,:)));
        for k = 1+n_cap:n_cap+n_per
            M = matMul2x2( M , T(:,:,k,:,:) );
            M = matMul2x2( M , R(:,:,k + 1,:) );
        end
        M = matPow2x2(M,n_ML);
        M_out = matMul2x2( M_out , M );
    end
    % Substrate layers
    for k = 1+n_cap+n_per:size(T,3)
        M_out = matMul2x2( M_out , T(:,:,k,:,:) );
        M_out = matMul2x2( M_out , R(:,:,k + 1,:,:) );
    end
    %%
    s = size(R);
    s(3) = [];
%     s([3 4]) = s([4 3]);
    M_out = reshape(M_out,s);
    %% Calculate reflected and transmitted field
    r = squeeze( M_out(2,1,:,:) ./ M_out(1,1,:,:) );
    t = squeeze( 1 ./ M_out(1,1,:,:) );
end

function A = matPow2x2(A,N)
    for i = 1:numel(A)/4
        A(:,:,i) = A(:,:,i)^N;
    end
end

function AB = matMul2x2(A,B)
    s = size(A);
    A = reshape(A,2,2,[]);
    B = reshape(B,2,2,[]);
    AB = zeros(size(A));
    
    AB(1,1,:) = A(1,1,:).*B(1,1,:) + A(1,2,:).*B(2,1,:);
    AB(1,2,:) = A(1,1,:).*B(1,2,:) + A(1,2,:).*B(2,2,:);
    AB(2,1,:) = A(2,1,:).*B(1,1,:) + A(2,2,:).*B(2,1,:);
    AB(2,2,:) = A(2,1,:).*B(1,2,:) + A(2,2,:).*B(2,2,:);
    
    AB = reshape(AB,s);
end