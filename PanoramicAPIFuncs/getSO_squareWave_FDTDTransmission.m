function scatteredOrders = getSO_squareWave_FDTDTransmission(simParams,illumParams,materialParams,orders)
    %% Loop through and run simulations
    scatteredOrders = zeros(length(orders),length(D),length(lambda),2,2);
    for i = 1:length(D)
        fprintf(['D ' num2str(i) ' out of ' num2str(length(D)) '\n'])
        tStart = tic;
        %% Load simulation with blank stack
        loadSetup('C:\Users\stuar\Documents\Panoramic\Scatterometry\test_stack_FBC_IdealReflector.sim');

        %% Set grid size
        setVariableValues('Nx',Nx);
        setVariableValues('dz',dz);
        %% Update the wavelength range and n/k
        setVariableValues('crTheta',theta);
        setVariableValues('wiz_wavelength',lambda);

        setVariableValues(['n' matStr],real(nk))
        setVariableValues(['k' matStr],-imag(nk))

        %% Set excitation and observation planes
        setVariableValues('zExc_FBC',z0)
        setVariableValues('zObs_FBC',0)
        %% Add absorber block across centered on domain with duty cycle D
        coordsX_tmp = (coordsX - mean(coordsX))*D(i) + mean(coordsX);
        coords = [coordsX_tmp coordsY z0 z0+t];
    %     z0 = z0+t;
        coordString = coords2String(coords);
        material = matStr;
        callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);

        if false
            %% Save setup
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','saveSetup','C:\Users\stuar\Documents\Panoramic\Scatterometry\test_IdealReflector.sim')
        end

        %% Run simulation
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','simulateMaskUsingTEMPESTpr2',1);

        %% Collect outputs
        [tableTE,~] = getEMOutputTable(outh,'Scattered TE Orders',true);
        [tableTM,header] = getEMOutputTable(outh,'Scattered TM Orders',true);

        %% Organize simulation outputs
        m = 0.5 + orders;
        tableTE(isnan(tableTE(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        tableTM(isnan(tableTM(:,strcmp(header,'n'))),strcmp(header,'n')) = 0;
        cgrp = unique(tableTE(:,strcmp(header,'cgrp')));

        [~,~,cols] = intersect({'cgrp','m'},header,'stable');
        dataCol = strcmp(header,'Data');

        %% Pull out the 0 order for each angle
        cgrp = [0 1];
        for iCGRP = 1:2
            for iO = 1:length(orders)
                %% TE Output
                rows = ismember(tableTE(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(iO,i,:,iCGRP,1) = tableTE(rows,dataCol);

                %% TM Output
                rows = ismember(tableTM(:,cols),[cgrp(iCGRP) 0.5+orders(iO)],'rows');
                scatteredOrders(iO,i,:,iCGRP,2) = tableTM(rows,dataCol);
            end
        end
        %%
        tEnd = toc(tStart);
        fprintf(['Finished in ' num2str(tEnd) ' s' '\n'])
    end
end