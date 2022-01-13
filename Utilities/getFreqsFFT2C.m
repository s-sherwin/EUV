function fx = getFreqsFFT2C(L,N)
    if mod(N,2) == 0
        fx = (-N/2:(N-2)/2)*L/N;
        %fx = ((-N + 2)/2:N/2)*L/N;
    else
        fx = ((-N + 1)/2:(N - 1)/2)*L/N;
    end
end