function [table,header] = getEMOutputTable(outh,outStr,ampPhase)
    if nargin < 3
        ampPhase = false;
    end
    %% Find output handles
    soh  = getDataSeriesHandlesForOutput(outh);
    
    useHandle = false(size(soh));
    for i = 1:length(soh)
        curSOH = soh(i);
        s = cell(callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getDataSeriesText',curSOH));
        
        useHandle(i) = strcmp(s{1},outStr);
    end
    
    %% Import the data for each usable handle
    useHandle = find(useHandle);
    curSOH = soh(useHandle(1));
    cHeaders = cell(callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getColumnHeaders',curSOH));
    n = length(cHeaders)+1;
     
    header = [cHeaders(:)', {'Data'}];
    table = zeros(0,n);
    
    if ~ampPhase
        for i = 1:length(useHandle)
            curSOH = soh(useHandle(i));
            y = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getDataAsColumns',curSOH);
            %% Create table form
            x = zeros(length(y)/n,n);
            for j = 1:n
                x(:,j) = y(j:n:end);
            end
            %% Append
            table = [table; x];
        end
    else
        for i = 1:2:length(useHandle)
            %% Amplitude
            curSOH = soh(useHandle(i));
            y = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getDataAsColumns',curSOH);
            %% Create table form
            x = zeros(length(y)/n,n);
            for j = 1:n
                x(:,j) = y(j:n:end);
            end
            %% Phase
            curSOH = soh(useHandle(i+1));
            y = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getDataAsColumns',curSOH);
            %% Multiply data columnn by phase
            x(:,end) = x(:,end).*exp(1i*y(n:n:end));
            %% Append
            table = [table; x];
        end
    end
end