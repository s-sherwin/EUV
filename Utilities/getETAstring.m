function out_strs = getETAstring(tellapsed,Ntot,Nellapsed,tUnit)
    if nargin < 4
        tUnit = 60;
    end
    switch tUnit
        case 1
            tStr = 'sec';
        case 60
            tStr = 'min';
        case 3600
            tStr = 'hrs';
        otherwise
            tStr = '';
    end
%     tellapsed = toc(tstart);

    tper = tellapsed/Nellapsed; % time per iteration
    
    ttot = tper*Ntot; % estimated total time assuming equal time per iteration
    trm = ttot - tellapsed; % estimated remaining time subtracting off ellapsed time
    
    % Output strings: percentage+time completed, and estimated ETA
    out_strs = {['Completed: ' num2str(Nellapsed/Ntot*100,'%10.1f') '%'],['Ellapsed: ' num2str(round(tellapsed/tUnit,1)) ' ' tStr],...
        ['Remaining: ' num2str(trm/tUnit,'%10.1f') ' ' tStr]};
end