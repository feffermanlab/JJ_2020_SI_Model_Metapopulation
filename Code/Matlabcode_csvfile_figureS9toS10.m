

function [t] = Matlabcode_csvfile_figureS9toS10(r,rd,mu,alpha,beta,MaxTime)  
% close all
%repeation for dealing with the uncertainty of initial values
rng(1,'twister'); %%fixed the random draw each time
s = rng;

if nargin == 0    
    
    r    = [0.7 0.7];                                               % natural growth rate for [genotype1  genotype2] %% this is 1 by 2 vector, 1 is robust and 2 is wild
    rd   = [0.07 0.07];                                             % growth rate of Is for [genotype1 genotype2]
    mu   = [0.25 0.25];                                             % natural death rate for [genotype1 genotype2]
    alpha   = [0.25 0.25*1.5];                                            % disease-introduce mortality rate for [genotype_robust genotype_wild]
    beta    = [0.0004 0.0004 0.0004 0.0004];                            % transmission rate among genotypes [11 12 21 22]
    MaxTime =60;  
    
end


N=5;  %number of patches 
m=4;   %number of equations


rin = 100+(220-100)*rand(N,1); %rin is initial value of susceptible for robust
win = 250+(320-250)*rand(N,1); %win is initial value for susceptible for wild type

rng(s);

X0 = [ rin(1);win(1);   2; 2];
kl = 3000;
for ii = 2:N
       %  [S0 robust        S0 wild      I0 robust   I0 robust ]% this is a 4 by 1 vector
    ini = [rin(ii);  win(ii);  0; 0]; %if the last two numbers equals to 0, that is the scenario when disease started from one local patch;
                                      %if the last two numbers equals to 2,
                                      %the same as the last two numbers in
                                      %X0, that represents the scenario
                                      %when disease started from all
                                      %patches.
    X0  = [X0; ini];
    
    kl  = [kl; 3000];  %carrying capacity of each patch 
   
end

%define migration matrix; 

mig=0:0.0125:0.05;

A1=[]; %S+I in R
A2=[]; %I in R+I in W
A3=[]; %S+I in W
A4=[]; %I in R
A5=[]; %I in wild

for k=1:length(mig)
mig_matrix0 = (ones(N)-eye(N))*0;      % the baseline matrix without migration 
%%each of migration matrix according to the topology in Fig. S8
mig_matrixa = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0;0 0 1 0 1;1 0 0 1 0]*mig(k);            
mig_matrixb = [0 1 1 0 0; 1 0 1 0 0; 1 1 0 1 1;0 0 1 0 1;0 0 1 1 0]*mig(k);         
mig_matrixc = [0 1 1 0 0; 1 0 1 1 0; 1 1 0 0 1;0 1 0 0 1;0 0 1 1 0]*mig(k);
mig_matrixd = [0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
mig_matrixe = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 0;0 1 1 0 0]*mig(k);
mig_matrixf = [0 1 0 0 1; 1 0 1 0 1; 0 1 0 1 1;0 0 1 0 1;1 1 1 1 0]*mig(k);
mig_matrixg = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
mig_matrixh = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 1;0 1 1 1 0]*mig(k);
mig_matrixi = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 1;1 0 1 0 1;0 1 1 1 0]*mig(k);
mig_matrixj = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 0;1 1 1 0 0]*mig(k);
mig_matrixk = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 1;1 1 1 1 0]*mig(k);
tspan = 0:1:MaxTime;                                % here 1 as time one step

mig_matrix = mig_matrixk; %the migration matrix, which can be modified based on which topological structure a-k is interested in

% pulse = 0;
% last_k = 0;
options =odeset('RelTol', 3e-8, 'AbsTol', 1e-15);
[t, pop]=ode45(@SI_eq, tspan, X0, options, r, rd, mu, alpha, beta, kl, N, mig_matrix0, mig_matrix, m);

% separate the results in different variables
S1=pop(:,1:m:m*N);
S2=pop(:,2:m:m*N);
I1=pop(:,3:m:m*N);
I2=pop(:,4:m:m*N);

S1_ave=mean(pop(:,1:m:m*N),2);
S2_ave=mean(pop(:,2:m:m*N),2);
I1_ave=mean(pop(:,3:m:m*N),2);
I2_ave=mean(pop(:,4:m:m*N),2);


A1=[A1 S1_ave+I1_ave];
A2=[A2 S2_ave+I2_ave];
A3=[A3 I1_ave+I2_ave];
A4=[A4 I1_ave];
A5=[A5 I2_ave];


end
dlmwrite('ave_SImodel_5patch_RT_matrixk_spread_1_5_20.csv',A1); % "matrixk" or "matrixa" based on adjusting #72, "spread" or "all patch" depending on the setup of the last two values of #34
dlmwrite('ave_SImodel_5patch_WT_matrixk_spread_1_5_20.csv',A2);
dlmwrite('ave_SImodel_5patch_infected_total_matrixk_spread_1_5_20.csv',A3);
dlmwrite('ave_SImodel_5patch_infected_RT_matrixk_spread_1_5_20.csv',A4);
dlmwrite('ave_SImodel_5patch_infected_WT_matrixk_spread_1_5_20.csv',A5);

end


function dPop = SI_eq(t, pop, r, rd, mu, alpha, beta, kl, N, mig_matrix0,mig_matrix, m)


    X = pop(1:m*N);
    dPop = zeros(m*N,1);

    % References:
    % X(1+m*j) = S1 in patch j+1
    % X(2+m*j) = I1 in patch j+1
    % X(3+m*j) = S2 in patch j+1
    % X(4+m*j) = I2 in patch j+1
    % X(5+m*j) = L1 in patch j+1
    ...
    
    for j = 0:N-1 
        Ni = X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j);    %population size at patch i: sum up all Ss and Is
        
        %define newborns for each genotype susceptibles
        %S1S1
        LS1S1 = check_zero(X(1+m*j)*X(1+m*j)/Ni*r(1)*(1-Ni/kl(j+1)));
%         
        %S1S2
        LS1S2 = check_zero(X(1+m*j)*X(2+m*j)/Ni*(r(1)+r(2))/2*(1-Ni/kl(j+1)));
        
        %S1I1
        LS1I1 = check_zero(X(1+m*j)*X(3+m*j)*2/Ni*(r(1)+rd(1))/2*(1-Ni/kl(j+1)));
        
        %S1I2
        LS1I2 = check_zero(X(1+m*j)*X(4+m*j)/Ni*(r(1)+rd(2))/2*(1-Ni/kl(j+1)));        
        
        %S2S2
        LS2S2 = check_zero(X(2+m*j)*X(2+m*j)/Ni*r(2)*(1-Ni/kl(j+1)));
        
        %S2I1
        LS2I1 = check_zero(X(2+m*j)*X(3+m*j)/Ni*(r(2)+rd(1))/2*(1-Ni/kl(j+1)));
          
        %S2I2
        LS2I2 = check_zero(X(2+m*j)*X(4+m*j)*2/Ni*(r(2)+rd(2))/2*(1-Ni/kl(j+1)));
        
        %I1I1
        LI1I1 = check_zero(X(3+m*j)*X(3+m*j)/Ni*rd(1)*(1-Ni/kl(j+1)));
        

        %I1I2
        LI1I2 = check_zero(X(3+m*j)*X(4+m*j)/Ni*(rd(1)+rd(2))/2*(1-Ni/kl(j+1)));

        %I2I2
        LI2I2 = check_zero(X(4+m*j)*X(4+m*j)/Ni*rd(2)*(1-Ni/kl(j+1)));
        

        %S1 in patch i
        tmpS1 = dot(mig_matrix0(:,j+1), X(1:m:end)) - X(1+m*j) * sum(mig_matrix0(:,j+1)) ;% migration term among patches with focal patch j
          
        dPop(1+j*m) =LS1S1 + LS1S2 + LS1I1 + LS1I2 + LS2I1 + LI1I1 + LI1I2 - beta(1)*X(1+m*j)*X(3+m*j)- beta(2)* X(1+m*j)*X(4+m*j)- mu(1)*X(1+m*j)+ tmpS1 ;
        
        %S2 in patch i
        tmpS2 = dot(mig_matrix(:,j+1), X(2:m:end)) - X(2+m*j) * sum(mig_matrix(:,j+1)) ;

        dPop(2+j*m)=LS1S2 + LS1I2 + LS2S2 + LS2I1 + LS2I2 + LI1I2 + LI2I2 - beta(3)*X(2+m*j)*X(3+m*j)- beta(4)* X(2+m*j)*X(4+m*j) - mu(2)*X(2+m*j) + tmpS2;
        
        %I1 in patch i
        tmpI1 = dot(mig_matrix0(:,j+1), X(3:m:end)) - X(3+m*j) * sum(mig_matrix0(:,j+1));
        
        dPop(3+j*m) = beta(1)*X(1+m*j)*X(3+m*j) + beta(2)* X(1+m*j)*X(4+m*j)-(alpha(1)+mu(1))*X(3+m*j)+ tmpI1;
        
        %I2 in patch i
        tmpI2 = dot(mig_matrix(:,j+1), X(4:m:end)) - X(4+m*j) * sum(mig_matrix(:,j+1)) ;
         
        dPop(4+j*m) = beta(3)*X(2+m*j)*X(3+m*j) + beta(4)* X(2+m*j)*X(4+m*j)-(alpha(2)+mu(2))*X(4+m*j) + tmpI2;%
       
    end
end

    function check0 = check_zero(pp)
        if pp < 0
            check0 = 0;
        else
        check0=pp;
        end
    end