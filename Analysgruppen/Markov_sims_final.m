%Markovsimulationer, körs TotSims gånger för populationsstorleken som
%specificeras i varje del-script. Beta, gamma och Tmax måste defineras till
%samma värde för samtliga simuleringar!

% Under simulationen kommer matriserna som skall innehålla varje simulation
% att bli större ifall det krävs. De initieras till en storlek och
% fördubblas sedan i storlek för varje gång deras dimensioner behöver
% utökas. Varje simulation lagrar S och I i matriser där varje rad är en
% ny simulation, kolonnerna utgörs därför av S(t) eller I(t).

% Koden sparar också datan i .mat-filer för att underlätta data analysen
% efteråt.

%%
clear all

%Parametarar:
N=1e+3;
I0=10;
beta=0.3;
gamma=0.2;
%Tidskala dagar:
Tmax=300;
n0=Tmax;

%Simulationsmatriser
TotSims=1e2;
I_each_sim=zeros(TotSims,n0);
S_each_sim=zeros(TotSims,n0);
TimeSteps=zeros(TotSims,n0);


tic
for i=1:TotSims

%Initierar variabler
    S=zeros(1,n0);
    S(1)=N-I0;
    I=zeros(1,n0);
    I(1)=I0;
    t=zeros(1,n0);
    j=1;
    N_resize=0;

while(I(j)>0 || t(j)<Tmax )
   a=beta*S(j)*I(j)/N;
   b=gamma*I(j);
   p1=a/(a+b);
   %Gillspie algoritmen: 
   u1=rand;
   u2=rand;
   
   
   if(u1<=p1)
       %infection:
       S(j+1)=S(j)-1;
       I(j+1)=I(j)+1;
       
   else
       %recovery:
       I(j+1)=I(j)-1;
       S(j+1)=S(j);
   end
   %Updating time vector with time for event:
   t(j+1)=t(j)-log(u2)/(a+b);
   j=j+1;
  
   
   %If the length of S,I,t is to small, increase the vectors size by *2. 
   if(j==length(S))
       temp_l=2*length(S);
       Sny=zeros(1,temp_l);
       Iny=zeros(1,temp_l);
       tny=zeros(1,temp_l);
       
       Sny(1:length(S))=S;
       Iny(1:length(S))=I;
       tny(1:length(S))=t;
       
       S=Sny;
       I=Iny;
       t=tny; 
       
       clear Sny Iny tny
       %Updates variable that keeps track of number of resizes required:
       N_resize=N_resize+1;
   end
   
end
    %If the number of coloumns had to increase during the simulation:
    if(j>size(S_each_sim,2))
       %Generating matrices with the new required dimensions:
       temp_l=length(S);
       S_each_ny=zeros(TotSims,j);
       I_each_ny=zeros(TotSims,j);
       TimeSteps_ny=zeros(TotSims,j);
         
       S_each_ny(1:TotSims,1:size(S_each_sim,2))=S_each_sim(:,:);
       I_each_ny(1:TotSims,1:size(S_each_sim,2))=I_each_sim(:,:);
       TimeSteps_ny(1:TotSims,1:size(S_each_sim,2))=TimeSteps(:,:);
       
       %Matrices are updated:
       S_each_sim=S_each_ny;
       I_each_sim=I_each_ny;
       TimeSteps=TimeSteps_ny;
       
       clear S_each_ny I_each_ny TimeSteps_ny
       
       %Append new elements: 
       S_each_sim(i,1:j)=S(1:j);
       I_each_sim(i,1:j)=I(1:j);
       TimeSteps(i,1:j)=t(1:j);
    else
        %Appends elements if the length of S,I or t are less than the colms
        %of S_each_sim /(I,t).
        temp_l=length(S);
        
        S_each_sim(i,1:j)=S(1:j);
        I_each_sim(i,1:j)=I(1:j);
        TimeSteps(i,1:j)=t(1:j); 
    end

end
 toc
%Save data:
save Markov_2_N_1e3.mat I_each_sim S_each_sim TimeSteps;

%%
clear all
disp('Simulation 1 complete: ')

%Parametrar
N=1e+4;
I0=10;
beta=0.3;
gamma=0.2;
%Tidskala dagar:
Tmax=300;
n0=Tmax;
TotSims=1e2;

%Simulationsmatriser
I_each_sim=zeros(TotSims,n0);
S_each_sim=zeros(TotSims,n0);
TimeSteps=zeros(TotSims,n0);


tic
for i=1:TotSims

%Initierar variabler
    S=zeros(1,n0);
    S(1)=N-I0;
    I=zeros(1,n0);
    I(1)=I0;
    t=zeros(1,n0);
    j=1;
    N_resize=0;

while(I(j)>0 || t(j)<Tmax )
   a=beta*S(j)*I(j)/N;
   b=gamma*I(j);
   p1=a/(a+b);
   %Gillspie algoritmen: 
   u1=rand;
   u2=rand;
   
   
   if(u1<=p1)
       %infection:
       S(j+1)=S(j)-1;
       I(j+1)=I(j)+1;
       
   else
       %recovery:
       I(j+1)=I(j)-1;
       S(j+1)=S(j);
   end
   %Updating time vector with time for event:
   t(j+1)=t(j)-log(u2)/(a+b);
   j=j+1;
  
   
   %If the length of S,I,t is to small, increase the vectors size by *2. 
   if(j==length(S))
       temp_l=2*length(S);
       Sny=zeros(1,temp_l);
       Iny=zeros(1,temp_l);
       tny=zeros(1,temp_l);
       
       Sny(1:length(S))=S;
       Iny(1:length(S))=I;
       tny(1:length(S))=t;
       
       S=Sny;
       I=Iny;
       t=tny;
       clear Sny Iny tny
       %Updates variable that keeps track of number of resizes required:
       N_resize=N_resize+1;
   end
   
end
    %If the number of coloumns had to increase during the simulation:
    if(j>size(S_each_sim,2))
       %Generating matrices with the new required dimensions:
       temp_l=length(S);
       S_each_ny=zeros(TotSims,j);
       I_each_ny=zeros(TotSims,j);
       TimeSteps_ny=zeros(TotSims,j);
         
       S_each_ny(1:TotSims,1:size(S_each_sim,2))=S_each_sim(:,:);
       I_each_ny(1:TotSims,1:size(S_each_sim,2))=I_each_sim(:,:);
       TimeSteps_ny(1:TotSims,1:size(S_each_sim,2))=TimeSteps(:,:);
       
       %Matrices are updated:
       S_each_sim=S_each_ny;
       I_each_sim=I_each_ny;
       TimeSteps=TimeSteps_ny; 
       clear S_each_ny I_each_ny TimeSteps_ny
       %Append new elements: 
       S_each_sim(i,1:j)=S(1:j);
       I_each_sim(i,1:j)=I(1:j);
       TimeSteps(i,1:j)=t(1:j);
    else
        %Appends elements if the length of S,I or t are less than the colms
        %of S_each_sim /(I,t).
        temp_l=length(S);
        
        S_each_sim(i,1:j)=S(1:j);
        I_each_sim(i,1:j)=I(1:j);
        TimeSteps(i,1:j)=t(1:j); 
    end

end
 toc
%Save data:
save Markov_2_N_1e4.mat I_each_sim S_each_sim TimeSteps;


%%
clear all
disp('Simulation 2 complete: ')

%Parametrar:
N=1e5;
I0=10;
beta=0.3;
gamma=0.2;
Tmax=200;
n0=Tmax;

%Simulationsmatriser
I_each_sim=zeros(TotSims,n0);
S_each_sim=zeros(TotSims,n0);
TimeSteps=zeros(TotSims,n0);

tic
for i=1:TotSims

%Initierar variabler
    S=zeros(1,n0);
    S(1)=N-I0;
    I=zeros(1,n0);
    I(1)=I0;
    t=zeros(1,n0);
    j=1;
    N_resize=0;

while(I(j)>0 || t(j)<Tmax )
   a=beta*S(j)*I(j)/N;
   b=gamma*I(j);
   p1=a/(a+b);
   %Gillspie algoritmen: 
   u1=rand;
   u2=rand;
   
   
   if(u1<=p1)
       %infection:
       S(j+1)=S(j)-1;
       I(j+1)=I(j)+1;
       
   else
       %recovery:
       I(j+1)=I(j)-1;
       S(j+1)=S(j);
   end
   %Updating time vector with time for event:
   t(j+1)=t(j)-log(u2)/(a+b);
   j=j+1;
  
   
   %If the length of S,I,t is to small, increase the vectors size by *2. 
   if(j==length(S))
       temp_l=2*length(S);
       Sny=zeros(1,temp_l);
       Iny=zeros(1,temp_l);
       tny=zeros(1,temp_l);
       
       Sny(1:length(S))=S;
       Iny(1:length(S))=I;
       tny(1:length(S))=t;
       
       S=Sny;
       I=Iny;
       t=tny; 
       
       clear Sny Iny tny
       %Updates variable that keeps track of number of resizes required:
       N_resize=N_resize+1;
   end
   
end
    %If the number of coloumns had to increase during the simulation:
    if(j>size(S_each_sim,2))
       %Generating matrices with the new required dimensions:
       temp_l=length(S);
       S_each_ny=zeros(TotSims,j);
       I_each_ny=zeros(TotSims,j);
       TimeSteps_ny=zeros(TotSims,j);
         
       S_each_ny(1:TotSims,1:size(S_each_sim,2))=S_each_sim(:,:);
       I_each_ny(1:TotSims,1:size(S_each_sim,2))=I_each_sim(:,:);
       TimeSteps_ny(1:TotSims,1:size(S_each_sim,2))=TimeSteps(:,:);
       
       %Matrices are updated:
       S_each_sim=S_each_ny;
       I_each_sim=I_each_ny;
       TimeSteps=TimeSteps_ny; 
       clear S_each_ny I_each_ny TimeSteps_ny
       %Append new elements: 
       S_each_sim(i,1:j)=S(1:j);
       I_each_sim(i,1:j)=I(1:j);
       TimeSteps(i,1:j)=t(1:j);
    else
        %Appends elements if the length of S,I or t are less than the colms
        %of S_each_sim /(I,t).
        temp_l=length(S);
        
        S_each_sim(i,1:j)=S(1:j);
        I_each_sim(i,1:j)=I(1:j);
        TimeSteps(i,1:j)=t(1:j); 
    end

end
 toc


save Markov_2_N_1e5.mat I_each_sim S_each_sim TimeSteps;

