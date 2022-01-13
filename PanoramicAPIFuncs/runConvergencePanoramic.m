function [xMin,struct] = runConvergencePanoramic(varStr,runStr,outStr,x0,alpha,errTol)
    %% Initialize vars
    struct = [];
    struct.x = x0;
    struct.t = [];
    struct.err = 0;
    %% Initial run
    tstart = tic;
    
    setVariableValues(varStr,x0); % Set variable
    outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient',runStr,1); % Run scattering simulation
    [table,~] = getEMOutputTable(outh,outStr,true); % Get output from simulation
    tableTrue = table;
    
    struct.t(end+1) = toc(tstart);
    %% Run until error tolerance is reached
    x = x0;
    curErr = 0;
    while curErr < errTol
        tstart = tic;
        
        x = x*alpha; % Update variable
        setVariableValues(varStr,x); % Set variable
        outh = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient',runStr,1); % Run scattering simulation
        [table,~] = getEMOutputTable(outh,outStr,true); % Get output from simulation
        if size(table,2) ~= size(tableTrue,2)
            break;
        end
        ii = ismember(tableTrue(:,1:end-1),table(:,1:end-1),'rows');
        r = reshape(table(:,end) - tableTrue(ii,end),[],1);
        curErr = max(abs(r));%r'*r;
        
        struct.err(end+1) = curErr;
        struct.x(end+1) = x;
        struct.t(end+1) = toc(tstart);
    end
    %% Set the best value of x: the 2nd to last one tried 
    if length(struct.x) > 1
        xMin = struct.x(end-1);
    else
        xMin = x0;
    end
    setVariableValues(varStr,xMin); % Set variable
end