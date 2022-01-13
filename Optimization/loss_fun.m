function cost = loss_fun(datCell,f,x,w,metric,alpha_avgpow,alpha_maxpow)
    if iscell(datCell)
        outputs = cell(size(datCell));
        [outputs{:}] = f(x);
        y = vec(cell2mat(vec(cellfun(@vec,datCell,'UniformOutput',false))));
        v = abs( cell2mat(vec(cellfun(@(y)reshape(y,[],size(x,2)),vec(outputs),'UniformOutput',false))));
    else
        y = datCell;
        if isnumeric(f)
            v = f;
        else
            v = f(x);
        end
    end
    
    if nargin < 4
        w = ones(size(y));
    end
    
    if nargin < 5
        metric = 'MSE';
    end
    
    nan_vals = isnan(y);
    y(nan_vals) = [];
    v(nan_vals) = [];
    w(nan_vals) = [];
    
    y = y.*w;
    v = v.*w;

    if strcmp(metric,'MSE')
        res = y - v;
    else % Scaling factor from least-squares projection
        res = zeros(size(v));
        for i = 1:size(v,2) % Multi-dimensional output: loop over columns
            res(:,i) = y - (v(:,i)\y)*v(:,i);
        end
    end
    
    cost = mean( abs(res).^2 ) ./ mean( abs(y).^2 );
%     cost = mean( abs(res).^2 );
    if nargin >= 6
        ymu = mean(reshape(y.^2,[],2,size(y,2)));
        vmu = mean(reshape(v.^2,[],2,size(v,2)));
        cost = cost + alpha_avgpow*reshape(nanmean(((ymu - vmu)./nanmean(ymu,[1 2])).^2,2),1,[]);
        
        if nargin >= 7
            w = ones(size(y,1),1);
            w(1:end/2) = 0;
            ymax = nanmax(reshape(w.*y.^2,[],2,size(y,2)));
            vmax = nanmax(reshape(w.*v.^2,[],2,size(v,2)));
            cost = cost + alpha_maxpow*reshape(nanmean(((ymax - vmax)./nanmean(ymu,[1 2])).^2,2),1,[]);
        end
    end
    
end