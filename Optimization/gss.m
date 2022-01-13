% Golden section search
function [x,fx] = gss(a,b,f,niter) % a: LB, b: UB, f: function
    gr = (sqrt(5)+1)/2; % Golden ratio
    fa = f(a);
    fb = f(b);
    if fa < fb
        x = a;
        fx = fa;
    else
        x = b;
        fx = fb;
    end
    c = b - (b - a) / gr; 
    d = a + (b - a) / gr;
    for i = 1:niter
        fc = f(c);
        fd = f(d);
        if fc < fd
            b = d;
            fb = fd;
            if fc < fa
                x = c;
                fx = fc;
            else
                x = a;
                fx = fa;
            end
        else
            a = c;
            fa = fc;
            if fd < fb
                x = d;
                fx = fd;
            else
                x = b;
                fx = fb;
            end
        end
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    end
end