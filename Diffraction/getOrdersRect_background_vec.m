function so = getOrdersRect_background_vec(x0,dx,T,T0,orders)
    %% Absorber block
    so = getOrdersRect_vec(x0,dx,T,orders);
    %% Background
    so = so + getOrdersRect_vec(x0+0.5,1-dx,T0,orders);
end