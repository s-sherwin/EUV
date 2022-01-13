function x = clip(x,LB,UB)
    x = max(x,LB);
    x = min(x,UB);
end