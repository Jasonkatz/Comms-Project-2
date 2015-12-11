function [ y ] = myintdump( x , T , low , high )
% intdump2( x , T , [min , max] )
% integrates and dumps over range min to max

    y = zeros(1,size(x,2)/T);
    for n = low:high
        y = x(n:T:end) + y;
    end

    y = y/(high-low);

end