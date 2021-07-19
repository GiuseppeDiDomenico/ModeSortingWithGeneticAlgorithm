function [result] = dominance(a, b)
flagA = 0;
flagB = 0;
[~,n] = size(a);

for i=1:n
    if (a(i) < b(i))
        flagA =  1;
    elseif (a(i) > b(i))
            flagB = 1;
    end 
end

if (flagA == flagB)
    result = 0;
elseif (flagA == 1)
    result = -1;
else
    result = +1;
end
end