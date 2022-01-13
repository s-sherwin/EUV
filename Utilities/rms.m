function y = rms(x,dim)
    if nargin == 1
        y = sqrt(mean(abs(x).^2));
    else
        y = sqrt(mean(abs(x).^2,dim));
    end
end