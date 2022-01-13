function res = fft2c(f)
    res = fftshift(fft2(ifftshift(f)));
%     res = ifftshift(fft2(fftshift(f)));
end
