function res = ifft2c(F)
    res = fftshift(ifft2(ifftshift(F)));
%     res = ifftshift(ifft2(fftshift(F)));
end