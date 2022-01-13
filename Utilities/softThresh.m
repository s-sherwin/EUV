function x = softThresh(y,lambda)
    x = zeros(size(y));
    ii = abs(y) > lambda;
    if ~isscalar(lambda)
        lambda = lambda(ii);
    end
    x(ii) = y(ii).*(abs(y(ii)) - lambda)./abs(y(ii)); % y(|y|-l)/|y|
    x(isnan(x)) = 0;
%     ii = y < -lambda;
%     x(ii) = y(ii) + lambda;
%     ii = y > lambda;
%     x(ii) = y(ii) - lambda;
end