function [f0f1,rho_nom,amu] = getScatFac(materials,lambda)
    dat = load('PeriodicTable.mat');
    
    [~,iA,iC] = intersect(materials,dat.elements,'stable');
    
    f0f1 = zeros(length(lambda),length(materials));
    amu = zeros(length(materials),1);
    rho_nom = zeros(length(materials),1);
    for i = 1:length(materials)
        ind = iC(i);
        
        amu(i) = dat.table.weight{ind};
        rho_nom(i) = dat.table.density{ind};
        
        ev = dat.table.asf{ind}(:,1);
        lam = ev2nm(ev);
        
        ff = dat.table.asf{ind}(:,2:3);
        
        f0f1(:,i) = interp1(lam,ff(:,1) + 1i*ff(:,2),lambda,'linear');
    end
%     f0f1 = dat.table.asf{inds(1)};
end