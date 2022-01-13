function so = getOrdersRect_vec(x0,dx,T,orders)
    orders = reshape(orders,[],1);
    so = zeros(length(orders),length(dx));
    %% Orders for centered rectangle of width dx (periodic with period 1)
    so(orders == 0,:) = dx;
    so(orders ~= 0,:) = sin(pi*orders(orders~=0).*dx)./(pi*orders(orders~=0));
    %% Scale by T (complex transmission)
    so = T.*so;
    %% Phase ramp for translation by x0
    so = so.*exp(-1i*2*pi*orders.*x0);
end