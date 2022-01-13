set(groot,'defaultAxesTickLabelInterpreter','latex');

filmModel = load('Reflectivityapp_workspace_fit131.mat');
filmModel = filmModel.model;
%%
lambda0 = 13.5;
theta0 = 6;
NA_mask = 0.33/4;
mag = 4;
NA_im = NA_mask*mag;
ux = linspace(-0.99,0.99,20); % Source points (x)
uy = ux'; % Source points (y)
[uxG,uyG] = meshgrid(ux,uy); % Source point grid
uxy = [uxG(:),uyG(:)];
uxy_init = uxy;
illum_in_pup = rssq(uxy_init,2)<=1;
uxy(~illum_in_pup,:) = []; % Remove points outside of pupil

uxy0 = [0,sind(theta0)/NA_mask];

lambda = lambda0;% + (-0.025:0.025:0.025)';

k0 = 2*pi./lambda;
kxy = k0.*reshape(uxy+uxy0,[1,size(uxy)])*NA_mask;
kz = sqrt(k0.^2 - sum(kxy.^2,3));
kxyz = [reshape(kxy,[],2),kz(:)];
% plot3(kxyz(:,1),kxyz(:,2),kxyz(:,3),'.')%,max(kxyz(:,1))*cos(0:0.01:2*pi),max(kxyz(:,2))*sin(0:0.01:2*pi),max(kxyz(:,3))*ones(size(0:0.01:2*pi)),'k')

illum_tbl = kxyz2illum(kxyz);

sigmaLambda = 0.1;%lambda0*0.01;
Slambda = exp(-0.5*((illum_tbl(:,1)-lambda0)/sigmaLambda).^2);

NA = 0.33/4;

k1 = 0.4; % k1
p = 2*k1*lambda0/NA_im; % pitch [wafer]
D = 0.5; % Duty cycle (same for xy)

z00 = 0;
zz = (-60:5:60)'; % Focus levels in FEM
DoF_target = 30; % Target for DoF
Wz = exp(-0.5*(zz/DoF_target).^2);

area_lattitude = 0.025;
exposure_lattitude = logspace(log10(1/1.1),log10(1.1),31) - 1; % Exposure levels offset from nominal in FEM
EL_target = 0.05; % Target for EL
We = exp(-0.5*(exposure_lattitude/EL_target).^2);
Wpw = Wz.*We; % Process window weights
I0metric = @(metric,thr) log(metric) - 0.01*log(thr); % Metric combining PW and TPT
st = strel('disk',0);

Nsc = 20; % Max scattered order index 
orders_x = -Nsc:Nsc; % Scattered orders (x)
orders_y = orders_x'; % Scattered orders (y)
theta_max = 60; % Max angle in calculation
sigmaF = 0.5*sind(theta_max)/(lambda0); % Apodization factor in scattering

%% Target pattern
D_buffer = 0;%.05; % Print buffer
D_buffer1 = 0;%.05; % Do not print buffer
P_print = double( max(abs(orders_x/length(orders_x)),abs(orders_y/length(orders_y))) <= (D-D_buffer)/2);
P_dontprint = double( max(abs(orders_x/length(orders_x)),abs(orders_y/length(orders_y))) >= (D+D_buffer1)/2);

%% Scattering
settings = [];
settings.model = filmModel;
% settings.illum_tbl = illum_tbl;
settings.orders_x = orders_x;
settings.orders_y = orders_y;

settings.Dx = D;
settings.Dy = D;
settings.theta_max = theta_max;
settings.sigmaF = sigmaF;


px = p;
py = p;
settings.px = px;
settings.py = py;

settings.model_str = 'DblSc';

settings.bandlim.theta_xy0 = [0, theta0]; % Mask-side CRA
settings.bandlim.NA = 0.5; % Bandlimit for scattering; does not need to equal image/mask bandlimit

%% Source
rad0y = lambda0/py/NA_im/2;
rad0x = lambda0/px/NA_im/2;
w = and(abs(uxy(:,1)) >= rad0x,abs(uxy(:,2)) >= rad0y);
kp_src = w > 0;

uxy_tmp = uxy(kp_src,:);

kp_src1 = vec(repmat(kp_src',length(lambda),1));

uxy_tmp = reshape(permute(repmat(uxy_tmp,1,1,length(lambda)),[3 1 2]),[],2);

illum_tbl_tmp = illum_tbl(kp_src1,:);
settings.illum_tbl = illum_tbl_tmp;
source0 = reshape(Slambda(kp_src1),1,1,[]); % Source, including only nonzero weights; including lambda weights


scatter_settings = [];

scatter_settings.settings = settings;
scatter_settings.illum_tbl = illum_tbl_tmp;

%% Near-field scattering
E = nearFieldDoubleScattering(settings);
%% Near-field for first source point
clf
subplot(121)
imagesc(abs(ifft2c(E(:,:,1))*numel(E(:,:,1))).^2)
colorbar
subplot(122)
imagesc(rad2deg(angle(ifft2c(E(:,:,1))*numel(E(:,:,1)))))
colorbar

%% Imaging settings
im_settings = [];
im_settings.NA_mask = NA_mask;
im_settings.NA_im = NA_im;
im_settings.mag_xy = [4 4];
im_settings.theta0 = theta0;
im_settings.zz = zz;
im_settings.z0 = z00;
im_settings.z0_fit = false;
im_settings.source = source0;

%% Exposure settings
exposure_settings = [];

exposure_settings.exposure_lattitude = exposure_lattitude;
exposure_settings.area_lattitude = area_lattitude;
exposure_settings.P_print = P_print;
exposure_settings.P_dontprint = P_dontprint;
exposure_settings.metric = 'Print no print';
exposure_settings.I0metric = I0metric;
exposure_settings.Wpw = Wpw;
exposure_settings.metric_PWerr = 'Mean err';
exposure_settings.st = st;

%% Partially coherent imaging 
[metric,outputs] = partialCohImMetric(scatter_settings,im_settings,exposure_settings);
%% FEM errormap
clf
imagesc(outputs.FEM)
colorbar

%% In-focus partially coherent image
clf
% imagesc(squeeze(outputs.I(:,1,:)))
subplot(131)
imagesc(outputs.I(:,:,1))
subplot(132)
imagesc(outputs.I(:,:,zz==0))
subplot(133)
imagesc(outputs.I(:,:,end))
