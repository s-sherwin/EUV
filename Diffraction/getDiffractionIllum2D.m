function illum_tbl_out = getDiffractionIllum2D(illum_tbl,orders_x,orders_y,px,py)
    kxyz = illum2kxyz(illum_tbl);
    
    k0 = 2*pi./illum_tbl(:,1);
    
    dkx = 2*pi*reshape(orders_x./px,1,[]);
%     dkx = repmat(dkx,
    dky = 2*pi*reshape(orders_y./py,1,1,[]);
    
    kx1 = kxyz(:,1) + dkx + 0*dky;
    ky1 = kxyz(:,2) + dky + 0*dkx;
    
    kx1 = clip(kx1,-repmat(k0,1,size(kx1,2)),repmat(k0,1,size(kx1,2)));
    ky1 = clip(ky1,-repmat(k0,1,size(ky1,2)),repmat(k0,1,size(kx1,2)));
%     can_prop = kx1.^2 + ky1.^2 < k0.^2;
    
    
    
    kz1 = real(sqrt(k0.^2 - kx1.^2 - ky1.^2));
    
    kxyz1 = cat(2,kx1(:),ky1(:),kz1(:));
    illum_tbl_out = kxyz2illum(kxyz1);
end