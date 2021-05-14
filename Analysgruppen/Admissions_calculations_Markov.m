%% Ladda in data:

clear all
clc

I0=10;
beta=0.078;
gamma=1/14;
%Tidskala dagar:
Tmax=25*7;
dt=1;
n0=int64(Tmax/dt);

TotSims=100;


%För att läsa in datan för olika befolkningsmängder:
% 1, 000 000 -befolkning:
%load Markov_1_N_1e6.mat
%N=1e+6;


% 100, 000 -befolkning
%load Markov_1_N_1e5.mat
%N=1e+5;

% 10, 000 -befolkning
%load Markov_1_N_1e4.mat
%N=1e+4;

% 1,000 -befolkning
%load Markov_1_N_1e3.mat
%N=1e+3;

load Markov.mat

%Interpolate data:
MaxTime=max(max(TimeSteps));
Time_interp=linspace(0,MaxTime,size(S_each_sim,2));
%I_interp=zeros(size(I_each_sim,1),size(I_each_sim,2));
S_interp=zeros(size(S_each_sim,1),size(S_each_sim,2));

%Interpolate each simulation as many times as the maxium number of time
%steps


for i=1:TotSims
    %Find first 0:element
    first_0=find(S_each_sim(i,:)==0, 1, 'first');
    
    if(first_0>0)
    %Fill time-steps with elements for each new time-step to interpolate:
    
    %t_final=linspace(TimeSteps(i,first_0-1),MaxTime+1,length(Time_interp)-length(TimeSteps(i,1:first_0-1))+1 );
    %TimeSteps(i,first_0:length(Time_interp))=t_final(2:end);
    %TimeSteps(i,first_0:length(Time_interp))=Time_interp(first_0:length(Time_interp));
    
    %Finaly do the linear interpolation:
    S_interp(i,:)=interp1(TimeSteps(i,1:first_0-1),S_each_sim(i,1:first_0-1),Time_interp,'linear','extrap');
    else
    S_interp(i,:)=interp1(TimeSteps(i,:),S_each_sim(i,:),Time_interp,'linear','extrap');
    end
    
end
    %Exchanges negative values with a zero :
    S_interp=max(S_interp,0);

%% Bilda medel-instans och mer 
clc

HospRate=0.05;
%Medelinstans:
%I_mean=mean(HospRate*I_interp);
S_mean=mean([zeros(size(S_each_sim,1),1),S_interp]);

% Standardavvikelse: 
% Sd_S=zeros(1,size(S_each_sim,1));
% Sd_I=zeros(1,size(I_each_sim,1));
% 
% for i=1:size(I_each_sim,1)
%     i-index, väljer rad
%     for j=1:n+1
%         j-index, väljer kolonn
%         Beräknar summan av (x_i - mu_i)^2 + ... +
%         Sd_S(1,i)=Sd_S(1,i)+(S_each_sim(i,j)-S_mean(1,j))^2;
%         Sd_I(1,i)=Sd_I(1,i)+(I_each_sim(i,j)-I_mean(1,j))^2;
%     end
%        
% end
%     Normaliserar standardavvikelsen:
%     Sd_S_norm=sqrt(Sd_S/double(n+1));
%     Sd_I_norm=sqrt(Sd_I/double(n+1));
    
%     %Identifiera instansen med maximal och minimal standardavviklese: 
%     [value_min,min_inst]=min(Sd_I_norm);
%     [value_max,max_inst]=max(Sd_I_norm);
    
   
    timestep=(Time_interp(end)-Time_interp(1))/length(Time_interp);
    Timesteps_OneWeek=int64(7/timestep);
    
    mean_adm_Week=zeros(1,25);
    %beräkna komulativa antalet inläggningar per vecka för 24 veckor med
    %startpunkt vecka 1 adm=0;
    for i=2:25
        if(Timesteps_OneWeek*25>length(S_mean))
        mean_adm_Week(i)=HospRate*abs(S_mean(end)-S_mean((i-1)*Timesteps_OneWeek));
        else
        mean_adm_Week(i)=HospRate*abs(S_mean(i*Timesteps_OneWeek)-S_mean((i-1)*Timesteps_OneWeek));
        end
    end
    
    mean_adm_komulativ=zeros(1,25);
    for i=1:24
        mean_adm_komulativ(i+1)=mean_adm_komulativ(i)+mean_adm_Week(i+1);
    end
    
    
    %Medel
    DeltaS_kom_mean=mean(mean_adm_komulativ);
    Sd_DeltaS_kom=0;
    %Standardavvikelse: 
    for i=1:length(mean_adm_komulativ)
        %Beräknar summan av (x_i - mu_i)^2 + ... +
        Sd_DeltaS_kom=Sd_DeltaS_kom+(mean_adm_komulativ(i)-DeltaS_kom_mean)^2;     
    end
       
    %Normaliserar standardavvikelsen:
    Sd_DeltaS_norm_kom=sqrt(Sd_DeltaS_kom/length(mean_adm_komulativ));
  
%% 95% konfidensintervall:
clc
set(0,'defaulttextinterpreter','latex');

%For a 95% confidence interval:
CI_t=1.645;

%confidence intervall is: mean +- CI_z*(sd/(sqrt(n))), n=antal
%observationer
n=length(mean_adm_komulativ);
I_CI_above= mean_adm_komulativ + CI_t*Sd_DeltaS_norm_kom/sqrt(n);
I_CI_below= mean_adm_komulativ - CI_t*Sd_DeltaS_norm_kom/sqrt(n);
 
%% Save for confidence interval plott:
% I_mean_N_1e3=I_mean;
% I_1_N_1e3= I_interp;
% CI_a_1e3= I_CI_above;
% CI_b_1e3= I_CI_below;
% t_int_N_1e3=Time_interp;
% 
% save 1e3.mat I_mean_N_1e3 CI_a_1e3 CI_b_1e3 t_int_N_1e3

% I_mean_N_1e4=I_mean;
% I_1_N_1e4= I_interp;
% CI_a_1e4= I_CI_above;
% CI_b_1e4= I_CI_below;
% t_int_N_1e4=Time_interp;

%save 1e4.mat I_mean_N_1e4 CI_a_1e4 CI_b_1e4 t_int_N_1e4

% I_mean_N_1e5=I_mean;
% I_1_N_1e5= I_interp;
% CI_a_1e5= I_CI_above;
% CI_b_1e5= I_CI_below;
% t_int_N_1e5=Time_interp;
% 
% save 1e5.mat I_mean_N_1e5 CI_a_1e5 CI_b_1e5 t_int_N_1e5


mean_adm_komulativ_Markov=mean_adm_komulativ;
CI_a_Markov= I_CI_above;
CI_b_Markov= I_CI_below;
t_int_Markov=Time_interp;

save Markov_plot.mat mean_adm_komulativ_Markov CI_a_Markov CI_b_Markov t_int_Markov

%% Figur ritare för figurer där spridningen jämförs för olika populationer, avnänds endast för I_mean datan. 
clear all 
clc

set(0,'defaulttextinterpreter','latex');
N1=1e3;
N2=1e4;
N3=1e5;
beta=0.3;
gamma=0.2;
R0=beta/gamma;

load 1e3.mat 
load 1e4.mat 
load 1e5.mat

lightblue = [170, 170, 255] / 255;
lightgreen = [170, 255, 170] / 255;
lightred = [255, 170, 170] / 255;

%1e3:
curve1a=CI_a_1e3/N1;
curve2a=CI_b_1e3/N1;
I_1e3=I_mean_N_1e3/N1;
t_1e3=t_int_N_1e3;

%1e4:
curve1b=CI_a_1e4/N2;
curve2b=CI_b_1e4/N2;
I_1e4=I_mean_N_1e4/N2;
t_1e4=t_int_N_1e4;

%1e5:
curve1c=CI_a_1e5/N3;
curve2c=CI_b_1e5/N3;
I_1e5=I_mean_N_1e5/N3;
t_1e5=t_int_N_1e5;


ttt1 = [t_1e3, fliplr(t_1e3)];
deviation_flip = [curve1a, fliplr(curve2a)];
fill(ttt1,deviation_flip, lightblue, 'LineStyle','none'); grid on; hold on;
loglog(t_1e3,I_1e3,'b'); hold on;

ttt2 = [t_1e4, fliplr(t_1e4)];
deviation_flip = [curve1b, fliplr(curve2b)];
fill(ttt2,deviation_flip, lightred, 'LineStyle','none');  hold on;
loglog(t_1e4,I_1e4,'r'); hold on;

ttt3 = [t_1e5, fliplr(t_1e5)];
deviation_flip = [curve1c, fliplr(curve2c)];
fill(ttt3,deviation_flip, lightgreen, 'LineStyle','none');  hold on;
loglog(t_1e5,I_1e5,'g'); hold on;
xlabel('Tidsenhet, dagar','FontSize',14)
ylabel('Andel infekterade, I(t)/N','FontSize',14)
title(['$\beta=$',num2str(beta),', $\gamma=$',num2str(gamma),', $\mathcal{R}_0=$',num2str(R0)],'Interpreter','latex','FontSize',14);  
legend('CI95: $N=10^3$','$\langle I(t_i) \rangle $','CI95: $N=10^4$','$\langle I(t_i) \rangle $','CI95: $N=10^5$','$\langle I(t_i) \rangle $','Interpreter','latex','FontSize',11)
axis([0 200 0 0.07]);


