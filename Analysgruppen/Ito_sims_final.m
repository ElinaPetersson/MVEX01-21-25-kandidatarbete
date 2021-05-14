%Simulationen av SDE:n k?rs TotSims g?nger. Det är inte lika komplex som
%Markov men lagrar sina värden pss som Markov gör i matriser, där varje rad
%är en ny simulation och kolonner motsvarar tidsteg.

clear all
clc

%Parametrar
TotSims=1e+6;
%Tidsskala: dagar
dt=1; 
%Stopptid:
T=200;
n=int64(T/dt);
N=1e6;
beta=0.2;
gamma=1/14; 
%Tidsskala: dagar
dt=1; 

sdt=sqrt(dt);
t=linspace(0,T,n+1);
I0=10;

I_each_sim=zeros(TotSims,n+1);
S_each_sim=zeros(TotSims,n+1);

for i=1:TotSims
    
    if(mod(i,1000)==0)
    disp(['Simulation nr:',num2str(i)])
    end
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
    %Varje rad motsvarar en simulation, medan varje kolonn ?r ett tidssteg
    I_each_sim(i,:)=I;
    S_each_sim(i,:)=S;

end

%% Ber?kna medelsimulation och standardavvikelsen f?r varje simulation till medlet.
clc
clf

%Medel trajectory:n
S_mean=mean(S_each_sim);
I_mean=mean(I_each_sim);

figure(1)
plot(t,I_mean,'r');
hold on
plot(t,S_mean,'bl');
legend('I(t), infected','S(t), succeptible')
xlabel('Tidsenhet, veckor')
ylabel('Antal')

Sd_S=zeros(1,TotSims);
Sd_I=zeros(1,TotSims);

for i=1:TotSims
    %i-index, v?ljer rad
    for j=1:n+1
        %j-index, v?ljer kolonn
        %Ber?knar summan av (x_i - mu_i)^2 + ... +
        Sd_S(1,i)=Sd_S(1,i)+abs((S_each_sim(i,j)-S_mean(1,j))^2);
        Sd_I(1,i)=Sd_I(1,i)+abs((I_each_sim(i,j)-I_mean(1,j))^2);
    end
    
    
end
    %Normaliserar standardavvikelsen:
    Sd_S_norm=sqrt(Sd_S/double(n+1));
    Sd_I_norm=sqrt(Sd_I/double(n+1));
    
    %Identifiera instansen med maximal och minimal standardavviklese: 
    [value_min,min_inst]=min(Sd_I_norm);
    [value_max,max_inst]=max(Sd_I_norm);
%% Runge-kutta F?r SIR-modellen:
clc
%Runge-kutta ordning 4/5 
tspan=[0,T];
x0=[N-I0,I0];
[t_RK,xsol]=ode45(@(t,x)SIR(t,x,beta,gamma,N),tspan,x0);
%plot(t,xsol(:,2));
I_RK=xsol(:,2);
S_RK=xsol(:,1);

%%   Plottar n?gra instanser av I
   clc
   
   figure(2)
   %Plottar instansen med maximal och minimal standard avvikelse samt medel instansen 
   plot(t,I_each_sim(min_inst,:),'r','linewidth',2)
   hold on 
   plot(t,I_each_sim(max_inst,:),'m','linewidth',2)
   hold on
   plot(t_RK,I_RK,'c','linewidth',2)
   hold on
   plot(t,I_mean,'b','linewidth',2)
   xlabel('Tidsenhet, veckor')
   ylabel('Antal')
   legend('\sigma_{min}','\sigma_{max}','SIR-numeriskt','medel')
   title(['Maximal och minimal standardavvikelse',' antalet simulationer=',num2str(TotSims)]);
   
   figure(3)
   for i=1:20
      sim=randi([1,TotSims]);
      hold on
      plot(t,I_each_sim(sim,:))
   end
   hold on 
   plot(t,I_mean,'b','linewidth',2)
   xlabel('Tidsenhet, veckor')
   ylabel('Antal')
   title('Några slumpm?ssiga instanser och medlet');
   
%% Save data:

save SDE_sims_1e6_N_1e6.mat I_each_sim S_each_sim;


%%
function dXdt = SIR(t,x,beta,gamma,N)
%Model ekvationerna
dXdt = [-beta*x(1)*x(2)/N ; beta*x(1)*x(2)/N-gamma*x(2)];
end