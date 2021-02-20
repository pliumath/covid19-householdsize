function Deltas = getDeltas(deltas,dDay,attack,l)
%Get a vector of parameter delta for every simulation day.

if length(deltas) ~= length(dDay) + 1
    
    error("Data mismatch.")
    
else

    TD = repmat(deltas(end),l,1);

    n = length(deltas);

    head = 1;

    for i = 1:(n-1)
    
        tail = dDay(i);
    
        TD(head:tail) = deltas(i);
        TD((tail+1):(tail+attack)) = (1:attack).*(deltas(i+1) - deltas(i))./attack + deltas(i);
    
        head = tail+attack+1;

    end

    Deltas = TD;

end
end

