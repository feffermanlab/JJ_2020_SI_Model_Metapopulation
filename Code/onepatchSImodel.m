%10/15/2018: To better understand the performance of the model, here JJ shrink the
%orginial 3 genotype code to 2 genotype based on the suggestions from Nina.
%Assumptions (simple version JJ finished on Oct.15,2018): N patches migration matrix indicating the connection and migration direction among those patches. Here, distance is not specifically be the reason to determine the migration rate. 
%Both genotypes have same demographic features in both patches

%%11/9/18-note of the meeting with Nina
%during this meeting, for conservation topic, we decided to do how the
%ratio of genotype a vs. b change with the ratio of Geno a mortality/Geno b
%mortality and patch number N. 

%the code would become 2genotype but can easily expand to N patches.
%Therefore, JJ changed the previous code (only one example N=3 with 3
%genotypes).

%JJ modified the initial values -- robust one is even smaller from 100 to
%150, while wild type is even larger 250-350. JJ does not think K matters
%too much, because as long as both types grow, they would reach K
%eventually.JJ also reduced the natural mortality in case migration does
%not function because host dies too fast.

%7-1-19: After modifying the newborn equations (previous equations lack 2
%for two different type combinations), we found that RT always win in one
%patch no matter what the initial values of RT and WT. In this case, even 
%including migration, I doubt at equlibrium, migration can prevent WT from 
%going to extinction. In other words, the Fig. 3, which studied how the ratio
%of equlibrium RT and WT changes with mortality ratio and WT migration may 
%end up as all infinite (or 0). This is because the current model is continuous, 
%which does not include the possibility of genetic drift (if one
%host genotype is in really low density, this genotype might go to
%extinction). Therefore, we should either write the model in discrete
%version or find a way to show the genetic drift in the continuous version.

%Here I used the cut off of the continuous version: if the newborn
%combinations have really low value, I assume that combinations go to
%extinction (to mimic the influence of genetic drift on mating).

%8-2-19: JJ added the situation when disease starts at one patch instead of
%all five all the same time.


function [t] = Newcode_2geno_5patch_evolution_10_26_19(r,rd,mu,alpha,beta,MaxTime)  
% close all
%repeation for dealing with the uncertainty of initial values
% for h=1
%the gradient of mortality ratio 11 numbers
mo=1:0.1:1.8;
A1=[]; %S+I in R
A2=[]; %I in R+I in W
A3=[]; %S+I in W
A4=[]; %I in R
A5=[];
for i=1:length(mo)
if nargin == 0                                                      %%if the number of input arguments is zero
    r    = [0.5 0.5];                                               % natural growth rate for [genotype1  genotype2] %% this is 1 by 2 vector
    rd   = [0.05 0.05];                                             % growth rate of Is for [genotype1 genotype2]
    mu   = [0.1 0.1];                                             % natural death rate for [genotype1 genotype2]
    alpha   = [0.1 0.1*mo(i)];                                            % disease-introduce mortality rate for [genotype1 genotype2]
    beta    = [0.0001 0.0001 0.0001 0.0001];                            % transmission rate among genotypes [11 12 21 22]
    MaxTime =100;                                                   % end of simulation (simulatin time)
end

% global pulse
% global last_k

N=1;  %number of patches 
m=4;   %number of equations

%initial conditions for N patches
%   S11 I11 S12 I12 L1 L4...; S21 I21 S22 I22 L1 L4...; S31 I31 S32 I32 L1
%   L4...
% X0 = [500;1;10;1;0;0;0;0;0;0;0;0;0;0;   10;1;500;1;0;0;0;0;0;0;0;0;0;0;   10;1;500;1;0;0;0;0;0;0;0;0;0;0];%let patch 1 has 500 S1 genotype 1; patch 2 has 500 S2 genotype 2.
% Kl = [1000;1000;1000];
% rng(h);
X0 = [180;  260;  1; 1];
kl = 3000;
% for ii = 1
%        %  [S0 geno1        S0 geno2      I0 geno1   I0 geno2  ]% this is a 4 by 1 vector
%     ini = [randi([150,200],1);  randi([250,300],1);  1; 1];% transpose(repelem(0,(m-4)))];
%     %ini = [300; 1;         300; 1;         300; 1];
%     X0  = [X0; ini];
%     
%     kl  = [kl; 3000];  %carrying capacity of each patch 
%    % Kl  = [Kl; 1000];  %carrying capacity of each patch
% end

%define migration matrix; there is no need to define migration from patch i
% according to Fig. 2 in the manuscript
% mig=0:0.00001:0.001;
% for k=1:length(mig)
% mig_matrix0 = (ones(N)-eye(N))*0;      % no migration of genotype 1 with small mortality
% mig_matrixa = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0;0 0 1 0 1;1 0 0 1 0]*mig(k);            
% mig_matrixb = [0 1 1 0 0; 1 0 1 0 0; 1 1 0 1 1;0 0 1 0 1;0 0 1 1 0]*mig(k);         
% mig_matrixc = [0 1 1 0 0; 1 0 1 1 0; 1 1 0 0 1;0 1 0 0 1;0 0 1 1 0]*mig(k);
% mig_matrixd = [0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
% mig_matrixe = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 0;0 1 1 0 0]*mig(k);
% mig_matrixf = [0 1 0 0 1; 1 0 1 0 1; 0 1 0 1 1;0 0 1 0 1;1 1 1 1 0]*mig(k);
% mig_matrixg = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
% mig_matrixh = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 1;0 1 1 1 0]*mig(k);
% mig_matrixi = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 1;1 0 1 0 1;0 1 1 1 0]*mig(k);
% mig_matrixj = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 0;1 1 1 0 0]*mig(k);
% mig_matrixk = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 1;1 1 1 1 0]*mig(k);
tspan = 0:1:MaxTime;                                % here 1 as time one step

mig_matrix = 0;
mig_matrix0 = 0;

% pulse = 0;
% last_k = 0;
options =odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, pop]=ode45(@SI_eq, tspan, X0, options, r, rd, mu, alpha, beta, kl, N, mig_matrix0, mig_matrix, m);
% X=[t,pop];
% dlmwrite('Geno2-5patch_SI_dynamics_trackingP-8-6-19.csv',X,'-append')
% separate the results in different variables
S1=pop(:,1:m:m*N);
S2=pop(:,2:m:m*N);
I1=pop(:,3:m:m*N);
I2=pop(:,4:m:m*N);

% X1 = [sum(S1)+sum(I1);sum(S2)+sum(I2)].';
%dlmwrite('2geno_1patch_S_R.csv',S1,'-append');
% dlmwrite('2geno_5patch_transient_matrixk-8-23-19_250generations_spreading_test1.csv',X1,'-append');
% plot everything
% subplot(2,3,1)
% h = plot(t,S1(:,1:N));
% legend(strcat('patch =',compose("%d",1:N)),'FontSize',6);
% legend boxoff
% %lgd.FontSize = 8;
% %axis([0, 100, 0, 50])
% xlabel 'Time';
% ylabel 'S1'
% 
% subplot(2,3,2) 
% h=plot(t,I1(:,1:N));
% xlabel 'Time';
% ylabel 'I1'
% 
% subplot(2,3,3) 
% h=plot(t,S1(:,1:N)+I1(:,1:N),'--');
% xlabel 'Time';
% ylabel 'S1+I1'
% 
% subplot(2,3,4)
% h=plot(t,S2(:,1:N)); 
% xlabel 'Time';
% ylabel 'S2'
% 
% subplot(2,3,5) 
% h=plot(t,I2(:,1:N));
% xlabel 'Time';
% ylabel 'I2'
% 
% subplot(2,3,6) 
% h=plot(t,S2(:,1:N)+I2(:,1:N),'--');
% xlabel 'Time';
% ylabel 'S2+I2'
A1=[A1 S1+I1];
A2=[A2 S2+I2];
A3=[A3 I1+I2];
A4=[A4 I1];
A5=[A5 I2];
end

dlmwrite('SImodel_onepatch_RT_total_10genelong12_10_19.csv',A1);
dlmwrite('SImodel_onepatch_WT_total_10genelong12_10_19.csv',A2);
dlmwrite('SImodel_onepatch_infected_total_10genelong12_10_19.csv',A3);
dlmwrite('SImodel_onepatch_infected_RT_10genelong12_10_19.csv',A4);
dlmwrite('SImodel_onepatch_infected_WT_10genelong12_10_19.csv',A5);


end
% end




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
    
    for j = 0
        Ni = X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j);    %population size at patch i: sum up all Ss and Is
        
        %define newborns for each genotype susceptibles
        %S1S1
        LS1S1 = check_zero(X(1+m*j)*X(1+m*j)/Ni*r(1)*(1-Ni/kl(j+1)));
%         if(X(1+m*j)*X(1+m*j)/(Ni*Ni)<0.05)
%             LS1S1=0.5*LS1S1;
%         end
%         
        %S1S2
        LS1S2 = check_zero(X(1+m*j)*X(2+m*j)/Ni*(r(1)+r(2))/2*(1-Ni/kl(j+1)));
        
%         if(X(1+m*j)*X(2+m*j)*2/(Ni*Ni)<0.05)
%             LS1S2=0.5*LS1S2;
%         end
        %S1I1
        LS1I1 = check_zero(2*X(1+m*j)*X(3+m*j)/Ni*(r(1)+rd(1))/2*(1-Ni/kl(j+1)));
%         if(X(1+m*j)*X(3+m*j)*2/(Ni*Ni)<0.05)
%             LS1I1=0.5*LS1I1;
%         end
        
        %S1I2
        LS1I2 = check_zero(X(1+m*j)*X(4+m*j)/Ni*(r(1)+rd(2))/2*(1-Ni/kl(j+1)));
        
%         if(X(1+m*j)*X(4+m*j)*2/(Ni*Ni)<0.05)
%             LS1I2=0.5*LS1I2;
%         end
        
        %S2S2
        LS2S2 = check_zero(X(2+m*j)*X(2+m*j)/Ni*r(2)*(1-Ni/kl(j+1)));
        
%         if(X(2+m*j)*X(2+m*j)/(Ni*Ni)<0.05)
%             LS2S2=0.5*LS2S2;
%         end
        
        %S2I1
        LS2I1 = check_zero(X(2+m*j)*X(3+m*j)/Ni*(r(2)+rd(1))/2*(1-Ni/kl(j+1)));
        
%         if(X(2+m*j)*X(3+m*j)*2/(Ni*Ni)<0.05)
%             LS2I1=0.5*LS2I1;
%         end
%         
        %S2I2
        LS2I2 = check_zero(X(2+m*j)*X(4+m*j)*2/Ni*(r(2)+rd(2))/2*(1-Ni/kl(j+1)));
        
%         if(X(2+m*j)*X(4+m*j)*2/(Ni*Ni)<0.05)
%             LS2I2=0.5*LS2I2;
%         end
        
        %I1I1
        LI1I1 = check_zero(X(3+m*j)*X(3+m*j)/Ni*rd(1)*(1-Ni/kl(j+1)));
        
%          if(X(3+m*j)*X(3+m*j)/(Ni*Ni)<0.05)
%             LI1I1=0.5*LI1I1;
%         end
        
        %I1I2
        LI1I2 = check_zero(X(3+m*j)*X(4+m*j)/Ni*(rd(1)+rd(2))/2*(1-Ni/kl(j+1)));
%         if(X(4+m*j)*X(3+m*j)*2/(Ni*Ni)<0.05)
%             LI1I2=0.5*LI1I2;
%         end
        
        %I2I2
        LI2I2 = check_zero(X(4+m*j)*X(4+m*j)/Ni*rd(2)*(1-Ni/kl(j+1)));
        
%         if(X(4+m*j)*X(4+m*j)/(Ni*Ni)<0.05)
%             LI2I2=0.5*LI2I2;
%         end
        
        
%         dlmwrite(fname,pulse,'-append')
        
%        global pulse
%        global last_k
       
       
%        tmpS1 = 0;
%        for k=1:100
%        if 500+(k-1)*500 < t &  t < 500+(k-1)*500+10 
%           
%             if k> last_k
%                 pulse = 0;
%             end 
%            if pulse == 0 
%                 pulse = 1;
%                 last_k = k;
%                 tmpS1 = 1000;
% %                 t
% %                 k
%        
% %            fname=sprintf('pulse_tracking_%d.csv',k);   
% %            dlmwrite(fname,pulse,'-append') ;   
% %            if size(csvread(fname))> 1
% %                pulse = 0;
% %                 
%            end
% %            tmpS2=10*pulse;  
%        end
%        end
        %S1 in patch i
        tmpS1 = dot(mig_matrix0(:,j+1), X(1:m:end)) - X(1+m*j) * sum(mig_matrix0(:,j+1)) ;% migration term among patches with focal patch j
       
        
%         tmp = LS1S1 + LS1S2 + LS1I1 + LS1I2 + LS2I1 + LI1I1 + LI1I2 - beta(1)*X(1+m*j)*X(3+m*j)- beta(2)* X(1+m*j)*X(4+m*j)- mu(1)*X(1+m*j)+ tmpS1 ; %=dS1/dt
%         dPop(1+j*m) = check_p(X(1+j*m),tmp);
%         
        dPop(1+j*m) =LS1S1 + LS1S2 + LS1I1 + LS1I2 + LS2I1 + LI1I1 + LI1I2 - beta(1)*X(1+m*j)*X(3+m*j)- beta(2)* X(1+m*j)*X(4+m*j)- mu(1)*X(1+m*j)+ tmpS1 ;
        
        %S2 in patch i
        tmpS2 = dot(mig_matrix(:,j+1), X(2:m:end)) - X(2+m*j) * sum(mig_matrix(:,j+1)) ;
%         tmpS2=0;
        
%         tmp = LS1S2 + LS1I2 + LS2S2 + LS2I1 + LS2I2 + LI1I2 + LI2I2 - beta(3)*X(2+m*j)*X(3+m*j)- beta(4)* X(2+m*j)*X(4+m*j) - mu(2)*X(2+m*j) + tmpS2; %=dS2/dt
%         dPop(2+j*m)=check_p(X(2+j*m),tmp);
%         
        dPop(2+j*m)=LS1S2 + LS1I2 + LS2S2 + LS2I1 + LS2I2 + LI1I2 + LI2I2 - beta(3)*X(2+m*j)*X(3+m*j)- beta(4)* X(2+m*j)*X(4+m*j) - mu(2)*X(2+m*j) + tmpS2;
        %I1 in patch i
        tmpI1 = dot(mig_matrix0(:,j+1), X(3:m:end)) - X(3+m*j) * sum(mig_matrix0(:,j+1));
%         tmpI1=0;
        
%         tmp = beta(1)*X(1+m*j)*X(3+m*j) + beta(2)* X(1+m*j)*X(4+m*j)-(alpha(1)+mu(1))*X(3+m*j)+ tmpI1;%=dI1/dt
%         dPop(3+j*m) = check_p(X(3+j*m),tmp);
        
        dPop(3+j*m) = beta(1)*X(1+m*j)*X(3+m*j) + beta(2)* X(1+m*j)*X(4+m*j)-(alpha(1)+mu(1))*X(3+m*j)+ tmpI1;
        %I2 in patch i
        tmpI2 = dot(mig_matrix(:,j+1), X(4:m:end)) - X(4+m*j) * sum(mig_matrix(:,j+1)) ;
         
%         tmpI2=0;
%         tmp = beta(3)*X(2+m*j)*X(3+m*j) + beta(4)* X(2+m*j)*X(4+m*j)-(alpha(2)+mu(2))*X(4+m*j) + tmpI2;%=dI2/dt
%         dPop(4+j*m) = check_p(X(4+j*m),tmp);
        
         dPop(4+j*m) = beta(3)*X(2+m*j)*X(3+m*j) + beta(4)* X(2+m*j)*X(4+m*j)-(alpha(2)+mu(2))*X(4+m*j) + tmpI2;%
       
%         %newborns Ls       
%         dPop(5+j*m) = X(1+m*j)*((X(1+m*j)-1)/(2*Ni))*r(1)*(1-Ni/Kl(j+1));   %babies from S1 & S1-L1
%         dPop(6+j*m) = X(3+m*j)*((X(3+m*j)-1)/(2*Ni))*r(2)*(1-Ni/Kl(j+1));   %babies from S2 & S2-L2
%         dPop(7+j*m) = X(1+m*j)*(X(3+m*j)/Ni)*((r(1)+r(2))/2)*(1/2)*(1-Ni/Kl(j+1)); %babies from S1 & S2 - L4
%         dPop(8+j*m) = X(1+m*j)*(X(2+m*j)/Ni)*((r(1)+rd(1))/2)*(1-Ni/Kl(j+1));%babies from S1 & I1 - L7
%         dPop(9+j*m) = X(1+m*j)*(X(4+m*j)/Ni)*((r(1)+rd(2))/2)*(1/2)*(1-Ni/Kl(j+1));%babies from S1 & I2 - L8
%         dPop(10+j*m) = X(3+m*j)*(X(2+m*j)/Ni)*((r(2)+rd(1))/2)*(1/2)*(1-Ni/Kl(j+1));%babies from S2 & I1 - L10
%         dPop(11+j*m) = X(3+m*j)*(X(4+m*j)/Ni)*((r(2)+rd(2))/2)*(1-Ni/Kl(j+1));%babies from S2 & I2 -L11
%         dPop(12+j*m) = X(2+m*j)*((X(2+m*j)-1)/(2*Ni))*rd(1)*(1-Ni/Kl(j+1)); %babies from I1 & I1 - L16
%         dPop(13+j*m) = X(4+m*j)*((X(4+m*j)-1)/(2*Ni))*rd(2)*(1-Ni/Kl(j+1)); %babies from I2 & I2 - L17
%         dPop(14+j*m) = X(2+m*j)*(X(4+m*j)/Ni)*((rd(1)+rd(2))/2)*(1/2)*(1-Ni/Kl(j+1)); %  babies from I1 & I2 - L19
    end
end

    function check0 = check_zero(pp)
        if pp < 0
            check0 = 0;
        else
        check0=pp;
        end
    end

%     function check_pos = check_p(p1,p2)
%         if p1+p2<0 | p1==0 
%            check_pos = -p1;
%         else
%         check_pos=p2;
%         end
%     end
