function [ims,mu_im] = mean_pad(ims00,W)
    mu_im = mean(ims00);
    ims = ims00.*W + mu_im.*(1 - W);
end    