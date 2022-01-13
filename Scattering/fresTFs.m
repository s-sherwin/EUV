function rt = fresTFs(model,illum_tbl,orders_x,orders_y,px,py,theta_max,theta_xy0,NA)
    if nargin < 8
        theta_xy0 = [0 0];
        NA = 1;
    end
    
    p = 1/(1/px + 1/py);
    illum_tbl2 = getDiffractionIllum2D(illum_tbl,orders_x,orders_y,px,py);
    kp = abs(illum_tbl2(:,2)) <= theta_max;
    kxyz = illum2kxyz(illum_tbl2);
    uxy = ((kxyz(:,1:2)./rssq(kxyz,2)) - sind(theta_xy0))./NA;
    kp = and(kp,rssq(uxy,2) <= 1);
    outputs = calcMLAbsEtchFresnel_model(model,illum_tbl2(kp,:));
    % outputs.r = cat(2,r_MLabs,r_MLetch,r_ML,r_abs,r_etch);
    % outputs.t = cat(2,t_abs,t_etch);
    % outputs.M = cat(4,M_ML,M_abs,M_etch,M_MLabs,M_MLetch);

%     r_MLetch = outputs.r(:,1);
%     r_MLabs = outputs.r(:,2);
%     r_ML = outputs.r(:,3);
%     r_abs = outputs.r(:,4);
%     r_etch = outputs.r(:,5);
%     t_abs = outputs.t(:,1);
%     t_etch = outputs.t(:,2);
    
    rt = zeros(size(illum_tbl2,1),7);
    rt(kp,:) = cat(2,outputs.r,outputs.t);
%     rt(kp,:) = cat(2,r_ML,r_abs,r_etch,t_abs,t_etch);
    rt = reshape(rt,size(illum_tbl,1),length(orders_y),length(orders_x),[]);
    rt = permute(rt,[2 3 1 4]);
%     rt = rt.*WLPF;
end