%Simulationen av SDE:n körs TotSims gånger (A).

clear all
clc

TotSims=1e2;
%Tidsskala: veckor
dt=1; 
%Stopptid 25veckor * 7 dagar/vecka:
T=25*7;
n=int64(T/dt);
I_each_sim=zeros(TotSims,n+1);
S_each_sim=zeros(TotSims,n+1);
N=1e4;
beta=0.093;
%Average 14 dagars tid 
gamma=1/14;
%Tidsskala: veckor
dt=1; 

sdt=sqrt(dt);
t=linspace(0,T,n+1);
I0=10;

for i=1:TotSims

    %Initierar variabler
    S=zeros(1,n);
    S(1)=N-I0;
    I=zeros(1,n);
    I(1)=I0;
    
        
    
        for j=1:n
            %Weiner processerna: 
            eta1=normrnd(0,1);
            eta2=normrnd(0,1);
            %Simulationen:
            S(j+1)= S(j)- beta*(S(j)*I(j))*dt/N -sqrt(beta*S(j)*I(j)/N)*sdt*eta1;
            I(j+1)= I(j) +  (beta*S(j)*I(j)/N -gamma*I(j))*dt + sqrt(beta*S(j)*I(j)/N)*sdt*eta1 -sqrt(gamma*I(j))*eta2*sdt;
            
                if(I(j+1)<0)
                    I(j+1)=0;
                end                 
        end

    %Lagrar I och S-kurvan vid varje tidsitteration som ett element i en matris.
    %Varje rad motsvarar en simulation, medan varje kolonn är ett tidssteg
    I_each_sim(i,:)=I;
    S_each_sim(i,:)=S;

end
%% Ladda in simuleringar (B)
clear all
clc

TotSims=1e+6;
%Tidsskala: veckor
dt=1; 
%Stopptid:
T=200;
n=int64(T/dt);
I_each_sim=zeros(TotSims,n+1);
S_each_sim=zeros(TotSims,n+1);
N=1e+6;
beta=0.2;
gamma=1/14;
%Tidsskala: veckor
dt=1;
I0=10;

sdt=sqrt(dt);
t=linspace(0,T,n+1);

load 'SDE_sims_1e6_N_1e6.mat'

%% Identifiera antalet fall då pandemin dog ut:
clc

%Söker igenom varje simulation:
No_outb=0;

for i=1:TotSims
    %Letar efter 0:or bland infected talen, för att kunna slänga bort dem och beräkna frekvensen 
    [NollExist,jj]=find(~I_each_sim(i,:));   
    %Om summan i NollExist>0 => det existerar 0:or  
    if(sum(NollExist)>0)
        %Utöka antalet:
        No_outb=No_outb+1;
        %Spara platsen:
        simNumber(i)=i;
    end

end

%Kan köras om fallen som pandemin dog skall exkluderas. 
%Remove_rows: elementen i vektorn innehåller raden i I_each_sim som dog.
 %[~,Remove_rows] = find(simNumber); %also returns vector v, which contains the nonzero elements of X.
 disp(['Antalet fall utan pandemi: ',num2str(No_outb)])
 
 %Nästa steg låt endast simulationerna som ej dör ut va kvar:
 %I_each_sim(Remove_rows,:)=[];
 %S_each_sim(Remove_rows,:)=[];
 
Prob_to_Die=No_outb/TotSims;
disp('********************************************************************************************************************')
disp('********************************************************************************************************************')
disp('********************************************************************************************************************')
disp('********************************************************************************************************************')
disp(['Sannolikheten för sjukdomen att dö ut utan pandemi: P=',num2str(Prob_to_Die)])
   
%% Beräknar medelinläggningar och sedan kumulativa medelinläggningar
clc
clf

HospRate=0.05;
%Medel trajectory:n
S_mean=mean(S_each_sim);
I_mean=mean([zeros(size(I_each_sim,1),1),HospRate*I_each_sim]);

    
    mean_adm_Week=zeros(1,25);
    %beräkna komulativa antalet inläggningar per vecka för 24 veckor med
    %startpunkt vecka 1 adm=0;
    for i=2:25
        mean_adm_Week(i)=HospRate*abs(S_mean((i)*7)-S_mean((i-1)*7));
    end
    
    mean_adm_komulativ=zeros(1,25);
    for i=1:24
        mean_adm_komulativ(i+1)=mean_adm_komulativ(i)+mean_adm_Week(i+1);
    end
    
    

  
%% Ta fram ett 95% konfidensintervall för medel-instansen:
clc
set(0,'defaulttextinterpreter','latex');

%For a 95% confidence interval:
CI_t=1.645;

%confidence intervall is: mean +- CI_t*(sd/(sqrt(n))), n=antal
%observationer
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
    
    %Konfidensintervallen:
    n=length(mean_adm_komulativ);
    I_CI_above= mean_adm_komulativ + CI_t*Sd_DeltaS_norm_kom/sqrt(n);
    I_CI_below= mean_adm_komulativ - CI_t*Sd_DeltaS_norm_kom/sqrt(n);



%%
mean_adm_komulativ_SDE=mean_adm_komulativ;
CI_a_SDE= I_CI_above;
CI_b_SDE= I_CI_below;
time_SDE=t;

save SDE_plot.mat mean_adm_komulativ_SDE CI_a_SDE CI_b_SDE time_SDE
   
