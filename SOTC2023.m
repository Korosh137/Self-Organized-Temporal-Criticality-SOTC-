% Korosh Mahmoodi 2016
% Mahmoodi, Korosh, Bruce J. West, and Paolo Grigolini. "Self-organizing complex networks: individual versus global rules."
% Frontiers in physiology 8 (2017): 478.


tic
clc ;
clear all ;
close all ;

TimeStep = 1e6 ; % number of trials; simulation length

% parameters of the Prisoner's Dilemma game
s = 0 ;   %    b < s + 1
b = 1.5 ; % b = 1 + 0.5 i.e. tempation to cheat is 0.5

K0 = 0 ; % initial value of control parameter of the agents
CC = 0.25 ; %  \khi parameter for updating control parameter of the agents (K values); connecting the two last payoffs (feedback) to the change of the agent's control parameter

g000 = 0.1  ; % parameter of the decision making rates

MM = 5 ; % number of agents on x-axis of the 2D lattice
NN = 7 ; % number of agents on y-axis of the 2D lattice
         % total number of agents = MM * NN 

Selfk1010 = zeros(TimeStep, 1) ; % k value of one of the agents
SelfkAve = zeros(TimeStep, 1) ;  % average of the k values

X = zeros(TimeStep, 1) ;   % mean field (considering states C as +1 and D as -1)
PiT = zeros(TimeStep, 1) ; % average benefit of the agents at each trial

SelfK1 = zeros(MM, NN) ;   % previous k values of the agents
Selfk2 = zeros(MM, NN) ;   % updated k values of the agents

Pi1 = zeros(MM, NN) ;  % previous payoffs values of the agents
Pi2 = zeros(MM, NN) ;  % updated payoffs of the agents

Decision = zeros(MM, NN) ;  % decisions of the agents ( C as 1 and D as 0 )

g = zeros(MM, NN) ; % number of neighbors of agents at state C (1)

% initially all the agents are defectors

for ti = 2 : TimeStep
   
    Alt2 = 0 ;
    Self2 = 0 ;
    for ii = 1 : MM
        for j1 = 1 : NN
            if Decision(ii, j1) == 1
                Alt2 = Alt2 + 1 ;
            end
            if Decision(ii, j1) == 0
                Self2 = Self2 + 1 ;
            end
        end
    end
    
    X(ti) = (Alt2 - Self2)/(MM*NN) ;
    
    % Counting 4 Neighbors    
    for i= 2 : (MM - 1)   
        for j= 2:(NN - 1)   
            g(i,j) = Decision(i-1, j) + Decision( i+1, j) + Decision( i, j-1) + Decision( i, j+1) ;
        end
    end
        
    for i= 2 : (MM - 1)
        g(i,1) = Decision( i-1, 1) + Decision( i+1, 1) + Decision( i, NN) + Decision( i, 2)  ;
    end
        
    for i= 2 : (MM - 1)
        g(i,NN) = Decision( i-1, NN) + Decision( i+1, NN) + Decision( i, (NN - 1)) + Decision( i, 1) ;
    end
        
    for j= 2 : (NN - 1)
        g(1,j) = Decision( MM, j) + Decision( 2, j) + Decision( 1, j-1) + Decision( 1, j+1) ;
    end
        
    for j= 2 : (NN - 1)
        g(MM,j) = Decision( (MM - 1), j) + Decision( 1, j) + Decision( MM, j-1) + Decision( MM, j+1) ;
    end
    
    g(1,1) = Decision( MM, 1) + Decision( 2, 1) + Decision( 1, NN) + Decision( 1, 2)  ;
    g(MM,1) = Decision( (MM - 1), 1) + Decision( 1, 1) + Decision( MM, NN) + Decision( MM, 2)  ;
    g(1,NN) = Decision( MM, NN) + Decision( 2, NN) + Decision( 1, (NN - 1)) + Decision( 1, 1) ;
    g(MM,NN) = Decision( (MM - 1), NN) + Decision( 1, NN) + Decision( MM, (NN - 1)) + Decision( MM, 1) ;
    
    %  Payoffs according to the Prisoner's Dilemma game
    
    for ii = 1 : MM
        for j2 = 1 : NN
            
            G = g(ii,j2) ;
            
            if Decision( ii, j2) == 1
                Pi2(ii,j2) = G * (1) + (4 - G) * (-s) ;
            else
                Pi2(ii,j2) = G * ( b) ;
            end
            PiT(ti ) = PiT(ti ) + Pi2(ii, j2) ;
        end
    end
        
    PiT(ti) = PiT(ti)/(1*MM*NN*1) ;  
    
    %  New Decisions 
    
    for i = 1 : MM
        for j = 1 : NN
            
            g0 = g(i,j) ;
            kc = SelfK1(i, j) ;
            
            BBBB = exp( kc * ( 1 - 0.5*g0   )) ;
            a1to0 = g000 * BBBB ;
            
            r = rand ;
            if Decision(i, j) == 1
                
                if (r <= a1to0 )
                    Decision(i, j) = 0 ;
                else
                    Decision(i, j) = 1 ;
                end
            else
                
                a0to1 = g000 / BBBB ;
                r = rand;
                
                if (Decision(i, j) == 0)
                    
                    if (r <= a0to1 )
                        Decision(i, j) = 1 ;
                    else
                        Decision(i, j) = 0 ;
                    end
                end
            end
                        
        end
    end
    
    SelfkAve(ti) = mean(mean(SelfK1)) ;  
    Selfk1010(ti) = SelfK1(2, 2) ;
    
        %  Updating the control parameters of the agents
    for ab = 1 : MM
        for ac = 1 : NN
            
            if  Pi1(ab, ac) ==  Pi2(ab, ac)
                Selfk2(ab, ac) =   1  *  SelfK1(ab, ac) ;
            else
                Selfk2(ab, ac) = SelfK1(ab, ac) + CC * (Pi2(ab, ac) - Pi1(ab, ac) )    / ((Pi2(ab, ac)) + (Pi1(ab, ac)) ) ;
            end
            
            if      Selfk2(ab, ac) <= 0     ||    g000 * exp(Selfk2(ab, ac))  > 1           ||     g000 / exp(Selfk2(ab, ac))  > 1
                Selfk2(ab, ac) = SelfK1(ab, ac) ;
            end
            
        end
    end
        
    for aa = 1 : MM
        for bb = 1 : NN
            SelfK1(aa, bb) = Selfk2(aa, bb) ;
            Pi1(aa, bb) = Pi2(aa, bb) ;
        end
    end
        
end  % end of trial

plot(SelfkAve)


toc