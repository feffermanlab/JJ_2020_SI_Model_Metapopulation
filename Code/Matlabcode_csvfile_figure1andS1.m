%%%this file can be modified to generate the csv datasets to generate Figures 1 and S1 


function [t] = Matlabcode_csvfile_figure1andS1(r,rd,mu,alpha,beta,MaxTime)  
% close all
%repeation for dealing with the uncertainty of initial values
% for h=1
%the gradient of mortality ratio 11 numbers
mo=1:0.5:2;
A1=[]; %S+I in R
A2=[]; %I in R+I in W
A3=[]; %S+I in W
A4=[]; %I in R
A5=[]; %I in W
for i=1:length(mo)
if nargin == 0                                                      %%if the number of input arguments is zero
    r    = [0.7 0.7];                                               % natural growth rate for [genotype_Robust  genotype_wild] %% this is 1 by 2 vector
    rd   = [0.07 0.07];                                             % growth rate of Is for [genotype_robust genotype_wild]
    mu   = [0.25 0.25];                                             % natural death rate for [genotype_robust genotype_wild]
    alpha   = [0.25 0.25*mo(i)];                                            % disease-introduce mortality rate for [genotype_robust genotype_wild]
    beta    = [0.0004 0.0004 0.0004 0.0004];                            % transmission rate among genotypes [11 12 21 22] "1" is robust, "2" is wild
    MaxTime =60;    
end


N=1;  %number of patches 
m=4;   %number of equations


X0 = [180;  250;  2; 2]; %arbitory initial values for S and I of both wild and robust  [S_robust, S_wild, I_robust, I_wild] 

kl = 3000; %carrying capacity

tspan = 0:1:MaxTime;                                % here 1 as time step

mig_matrix = 0;
mig_matrix0 = 0;


options =odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%%%call the function SI_eq to run the ode
[t, pop]=ode45(@SI_eq, tspan, X0, options, r, rd, mu, alpha, beta, kl, N, mig_matrix0, mig_matrix, m);

% separate the results in different variables
S1=pop(:,1:m:m*N); %for Susceptibles in Wild
S2=pop(:,2:m:m*N); %for susceptibles in robust
I1=pop(:,3:m:m*N); %for infected in wild
I2=pop(:,4:m:m*N); %for infected in robust


A1=[A1 S1+I1];
A2=[A2 S2+I2];
A3=[A3 I1+I2];
A4=[A4 I1];
A5=[A5 I2];
end

dlmwrite('SImodel_onepatch_RT_total_15gene.csv',A1);
dlmwrite('SImodel_onepatch_WT_total_15gene.csv',A2);
dlmwrite('SImodel_onepatch_infected_total_15gene.csv',A3);
dlmwrite('SImodel_onepatch_infected_RT_15gene.csv',A4);
dlmwrite('SImodel_onepatch_infected_WT_15gene.csv',A5);


end


%SI_eq function of the whole SI ode model structure
function dPop = SI_eq(t, pop, r, rd, mu, alpha, beta, kl, N, mig_matrix0,mig_matrix, m)
    %"1" is robust, "2" is for wild type

    X = pop(1:m*N);
    dPop = zeros(m*N,1);

    % References:
    % X(1+m*j) = S1 in patch j+1
    % X(2+m*j) = I1 in patch j+1
    % X(3+m*j) = S2 in patch j+1
    % X(4+m*j) = I2 in patch j+1
    % X(5+m*j) = L1 in patch j+1
    ...
    
    for j = 0
        Ni = X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j);    %population size at patch i: sum up all Ss and Is
        
        %define newborns for each genotype susceptibles
        %S1S1
        LS1S1 = check_zero(X(1+m*j)*X(1+m*j)/Ni*r(1)*(1-Ni/kl(j+1)));
                 
        %S1S2
        LS1S2 = check_zero(X(1+m*j)*X(2+m*j)/Ni*(r(1)+r(2))/2*(1-Ni/kl(j+1)));
        
        %S1I1
        LS1I1 = check_zero(2*X(1+m*j)*X(3+m*j)/Ni*(r(1)+rd(1))/2*(1-Ni/kl(j+1)));
        
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
%         
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
