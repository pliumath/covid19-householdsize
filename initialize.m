function B = initialize(A,E,I,R)
%Initialize individuals' states: 1-susceptible, 2-exposed, 3-infected,
%4-removed.

if size(A,1) < E + I + R
    
    error("More suspectible agents are needed.")
    
else
    
    T = A;
    T(:,2) = 1;
    
    x = randsample((1:size(A,1))',E);
    T(x,2) = 2;
    
    y = randsample(setdiff((1:size(A,1))',x),I);
    T(y,2) = 3;
    
    z = randsample(setdiff((1:size(A,1))',[x;y]),R);
    T(z,2) = 4;
    
    B = T;
    
end

end

