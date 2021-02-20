function [A,LH] = generateIndividuals(SH,HD,Q)
%Generate individuals and their corresponding household information.

if nargin < 3
    
    L = gendist(HD',SH,1);
    
    T = [];
    n = 1;
    
    for i = 1:SH
        
        m = n+L(i)-1;
        t = [repmat(i,L(i),1),zeros(L(i),1),(1:L(i))',repmat(L(i),L(i),1)];
        T = [T;t];
        n = m+1;
        
    end
    
    A = T;
    LH = count(L);
    
else
    
    L = gendist(HD',SH,1);
    
    T = [];
    n = 1;
    
    for i = 1:SH
        
        QP = gendist(Q,L(i),1);
        
        m = n+L(i)-1;
        t = [repmat(i,L(i),1),zeros(L(i),1),(1:L(i))',repmat(L(i),L(i),1),QP,zeros(L(i),1),randsample(2:4,L(i),1)'];
        T = [T;t];
        n = m+1;
        
    end
    
    A = T;
    LH = count(L);
    
end

end

function C = count(T)

C = [sum(T==1),sum(T==2),sum(T==3),sum(T==4),sum(T==5),sum(T==6),sum(T==7),sum(T==8)];

end

