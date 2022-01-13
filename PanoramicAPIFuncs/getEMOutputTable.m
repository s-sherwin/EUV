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
        
        useHandle(i) = strncmp(s{1},outStr,length(outStr));
    end
    
    %% Import the data for each usable handle
    useHandle = find(useHandle);
    curSOH = soh(useHandle(1));
    cHeaders = cell(callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getColumnHeaders',curSOH));
%     n = length(cHeaders)+1;
     
%     header = [cHeaders(:)', {'Data'}];
    headerCell = cell(length(length(useHandle)),1);
    tableCell = cell(length(length(useHandle)),1);
%     table = zeros(0,n);
    
    if ~ampPhase
        for i = 1:length(useHandle)
            curSOH = soh(useHandle(i));
            cHeaders = cell(callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getColumnHeaders',curSOH));
            n = length(cHeaders)+1;
            y = callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getDataAsColumns',curSOH);
            %% Create table form
            x = zeros(length(y)/n,n);
            for j = 1:n
                x(:,j) = y(j:n:end);
            end
            %% Append
            headerCell{i} = [cHeaders(:)', {'Data'}];
            tableCell{i} = x;
%             table = [table; x];
            
        end
    else
        for i = 1:2:length(useHandle)
            %% Amplitude
            curSOH = soh(useHandle(i));
            cHeaders = cell(callJavaFunction('panoramictech.v700.OpenAPI.CAPIClient','getColumnHeaders',curSOH));
            n = length(cHeaders)+1;
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
            headerCell{i} = [cHeaders(:)', {'Data'}];
            tableCell{i} = x;
%             table = [table; x];
        end
    end
    
    %% Stitch together cells into one table (because some may not have the same columns)
    headerCell(cellfun(@isempty,headerCell)) = [];
    tableCell(cellfun(@isempty,tableCell)) = [];
    %% header: union of all headers
    header = headerCell{1};
    for i = 1:length(headerCell)
        header = union(header,headerCell{i},'stable');
    end
    %% Now fill in missing columns
    tableCell_full = cell(size(tableCell));
    table = zeros(0,length(header));
    for i = 1:length(tableCell)
        %%
        tableCell_full{i} = NaN(size(tableCell{i},1),length(header));
        [~,iA,iB] = intersect(header,headerCell{i},'stable');
        tableCell_full{i}(:,iA) = tableCell{i}(:,iB);
        %%
        table = [table; tableCell_full{i}];
    end
end