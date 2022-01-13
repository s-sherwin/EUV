function kxyz = illum2kxyz(illum_tbl)
    k0 = 2*pi./illum_tbl(:,1);
    theta = illum_tbl(:,2);
    phi = illum_tbl(:,3);
    
    kz = k0.*cosd(theta);
    kx = k0.*sind(theta).*cosd(phi);
    ky = k0.*sind(theta).*sind(phi);
    
    kxyz = cat(2,kx,ky,kz);
end