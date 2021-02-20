%% A Discrete-time individual-based household-stratified SEIR model
% Code for the second scenario described in Modelling the impact of
% household size distribution on the transmission dynamics of COVID-19

%% Preparation

clear,clc

load('Distributions.mat') %load household size distributions

%rng(1); %fix seed

%% Initialization

SH = 5000; %initial susceptible households

Region = 1; %region: 1-VCH, 2-FH.
HD = HD(:,Region);

Q = [0.45,0,0.55]; %isolation preference: no isolation, home isolation, separated isoaltion

[A,LH] = generateIndividuals(SH,HD,Q);

N = size(A,1); %initial number of individuals
E(1) = 50; %initial number of the exposed
I(1) = 20; %initial number of the infected
R(1) = 0; %initial number of the removed
S(1) = N - E(1) - I(1) - R(1); %initial number of the susceptible

A = initialize(A,E,I,R);

%% Parameters

beta1 = 0.011; %community transmission probability
beta2 = 0.09; %household transmission probability
u = 0.15; %probability E -> I
v = 0.071; %probability I -> R

rho = 20; %number of contacts

deltas = [0.675,0.675]; %social distancing parameter: 0-never staying home, 1-always staying home
dDay =       [0]; %the days social distancing changes
attack = 7; %transforamtion time from a value of delta to the next one

qs = 0; %initial isolation status
qDay = 50; %the day isolation starts

l = 200; %length of simulation

isolation_delay = false; %scenario options
separated_isolation_offer = false;

%% Simulation

S(2:l,1) = 0; %preallocate variables
E(2:l,1) = 0;
I(2:l,1) = 0;
R(2:l,1) = 0;

C(1) = 0; %number of community transmissions
H(1) = 0; %number of household transmissions

RH = zeros(8,1); %recovered househol counts

Records = {A}; %record individuals' states

Deltas = getDeltas(deltas,dDay,attack,l); %get everyday parameter delta

for i = 1:l
    
    i
    
    delta = Deltas(i);
    
    if i > qDay
        qs = 1;
    end
    
    [S(i+1,1),E(i+1,1),I(i+1,1),R(i+1,1),Records{i+1,1},C,H,LH,RH] = update(Records{i,1},C,H,HD,LH,RH,beta1,beta2,u,v,delta,rho,Q,qs,isolation_delay,separated_isolation_offer);
    
end

%% Results and Figures

[InfCurv,IncPrd,InfPrd] = getResults(Records,LH,2);

makeFigures(I,C,H,InfCurv,HD,l)

%% Equations

function [S,E,I,R,B,C,H,LH,RH] = update(A,C,H,HD,LH,RH,beta1,beta2,u,v,delta,rho,Q,qs,isolation_delay,separated_isolation_offer)

T = A;
N = size(A,1);

TC = 0;
TH = 0;

for i = 1:N %check on every individual
    
    switch T(i,2) %determine the individual's state
        
        case 1 %if the individual is susceptible
            
            M = findHouseholdMembers(i,A);
            
            m = sum(M(:,2) == 3); %number of infected household members
            n = sum(A(:,2) == 3); %number of infected individuals
            
            if isolation_delay
                
                if qs == 1
                    Ms = M(M(:,2) == 3,5:7);
                    msq = sum(Ms==3&Ms(:,2)>=Ms(:,3));%number of household membesr separately isolated
                    mnq = sum(Ms(:,1)==1)+sum(Ms(:,1)~=1&Ms(:,2)<Ms(:,3));%number of household members not isolated
                    
                    As = A(A(:,2) == 3,5:7);
                    nnq = sum(As(:,1)==1)+sum(As(:,1)~=1&As(:,2)<As(:,3));%number of infected not isoalted
                    
                    c = (1-delta)*beta1*rho*((nnq-mnq)/N);
                    h = delta*beta2*(m-msq);
                else
                    c = (1-delta)*beta1*rho*((n-m)/N);
                    h = delta*beta2*m;
                end
                
            else
                
                if qs == 1
                    Ms = M(M(:,2) == 3,5);
                    msq = sum(Ms==3);%number of household members separately isolated
                    mnq = sum(Ms==1);%number of household members not isolated
                    
                    As = A(A(:,2) == 3,5);
                    nnq = sum(As==1);%number of infected not isolated
                    
                    c = (1-delta)*beta1*rho*((nnq-mnq)/N);
                    h = delta*beta2*(m-msq);
                else
                    c = (1-delta)*beta1*rho*((n-m)/N);
                    h = delta*beta2*m;
                end
                
            end
            
            p = c + h;  %community transmission and household transmission
            
            r = rand;
            
            if r < p
                
                T(i,2) = 2;
                
                if r < c %compute the number of community transmissions
                    
                    TC = TC+1;
                    
                else %compute the number of household transmissions
                    
                    TH = TH+1;
                    
                end
                
                
            end
            
        case 2  %if the individual is exposed
            
            if rand < u
                
                T(i,2) = 3;
                
            end
            
        case 3 %if the individual is infected
            
            T(i,6) = T(i,6)+1;
            
            if rand < v
                
                T(i,2) = 4;
                
            end
            
        case 4 %in the individual is recovered
            
            M = findHouseholdMembers(i,A);
            nhou = M(1,1);
            nmem = M(1,4);
            m = sum(M(:,2) == 4);
            
            if m == nmem %if all household members are recovered
                
                RH(nmem) = RH(nmem)+1;
                T(T(:,1) == nhou,2) = 0;
                
            end
            
    end
    
end

if separated_isolation_offer
    
    T(T(:,6)>T(:,7)+2 & T(:,5)==2 & T(:,4)>=5,5) = 3;
    
end

B = T; %updating states and results

S = sum(B(:,2)==1);
E = sum(B(:,2)==2);
I = sum(B(:,2)==3);
R = sum(B(:,2)>=4)+(1:8)*RH;

C = [C;TC];
H = [H;TH];

while R > 0.1*N
    
    NN = gendist(HD',1,1);
    NH = B(end,1);
    
    LH(NN) = LH(NN)+1;
    
    QP = gendist(Q,NN,1);
    
    for i = 1:NN
        B(end+1,:) = [NH+1,1,i,NN,QP(i),0,randsample(2:4,1,1)'];
    end
    
    N = size(B,1);
    
end

end