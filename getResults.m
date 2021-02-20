function [InfCurv,IncPrd,InfPrd] = getResults(Records,LH,Scn)
%Get results from the simulation.

if Scn == 1
    
    BC = Records;
    t = size(BC,1); %number of days simulated
    n = size(BC{end},1); %number of individuals
    
    B = zeros(n,4,t);
    for i = 1:t
        
        B(1:size(BC{i},1),:,i) = BC{i};
        
    end
    
    T = B(:,2,:);
    T = reshape(T,n,t);
    
    R = zeros(n,4);
    
    for i = 1:n
        
        for j = 1:4
            
            R(i,j) = sum(T(i,:)==j);
            
        end
        
    end
    
    IncPrd = mean(R(R(:,2)~=0,2));
    InfPrd = mean(R(R(:,3)~=0,3));
    
    PP = LH.*(1:8).*(500000/sum(LH));
    
    ICN = ones(t,8);
    ICD = ones(t,8);
    IC = ones(t,8);
    
    for i = 1:t
        
        BB = B(:,:,i);
        
        for j = 1:8
            
            Bj = BB(BB(:,4) == j,:);
            
            ICD(i,j) = sum(Bj(:,2) >= 3)+sum(Bj(:,2) == 0); %number of infected individuals from households of size j up to day i
            ICN(i,j) = PP(j) - ICD(i,j); %number of uninfected individuals from households of size j up to day i
            
        end
        
    end
    
    ICD2 = [zeros(1,8);ICD(1:(end-1),:)];
    ICD = ICD-ICD2;
    
    
    ICR = 1-ICD./ICN; %Kaplan-Meier estimator
    for i = 1:t
        
        IC(i,:) = prod(ICR(1:i,:));
        
    end
    IC(1,:) = 1;
    
    InfCurv = IC;
    
else
    
    BC = Records;
    t = size(BC,1); %number of days simulated
    n = size(BC{end},1); %number of individuals
    
    B = zeros(n,7,t);
    for i = 1:t
        
        B(1:size(BC{i},1),:,i) = BC{i};
        
    end
    
    T = B(:,2,:);
    T = reshape(T,n,t);
    
    R = zeros(n,4);
    
    for i = 1:n
        
        for j = 1:4
            
            R(i,j) = sum(T(i,:)==j);
            
        end
        
    end
    
    IncPrd = mean(R(R(:,2)~=0,2));
    InfPrd = mean(R(R(:,3)~=0,3));
    
    PP = LH.*(1:8).*(500000/sum(LH));
    
    ICN = ones(t,8);
    ICD = ones(t,8);
    IC = ones(t,8);
    
    for i = 1:t
        
        BB = B(:,:,i);
        
        for j = 1:8
            
            Bj = BB(BB(:,4) == j,:);
            
            ICD(i,j) = sum(Bj(:,2) >= 3)+sum(Bj(:,2) == 0); %number of infected individuals from households of size j up to day i
            ICN(i,j) = PP(j) - ICD(i,j); %number of uninfected individuals from households of size j up to day i
            
        end
        
    end
    
    ICD2 = [zeros(1,8);ICD(1:(end-1),:)];
    ICD = ICD-ICD2;
    
    
    ICR = 1-ICD./ICN; %Kaplan-Meier estimator
    for i = 1:t
        
        IC(i,:) = prod(ICR(1:i,:));
        
    end
    IC(1,:) = 1;
    
    InfCurv = IC;
end

end

