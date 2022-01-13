function y = rssq(x,dim)
    if nargin == 1
        y = sqrt(sum(abs(x).^2));
    else
        y = sqrt(sum(abs(x).^2,dim));
    end
end