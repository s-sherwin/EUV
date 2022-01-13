function illum_tbl = kxyz2illum(kxyz)
    kx = kxyz(:,1);
    ky = kxyz(:,2);
    kz = kxyz(:,3);
    k0 = sqrt(sum(kxyz.^2,2));
    lambda = 2*pi./k0;
    
    theta = real(acosd(kz./k0));
    phi = atan2d(ky,kx);
    
    illum_tbl = cat(2,lambda,theta,phi);
end