%Part 2  

close all
clear all
clc 
%Atencao, demora cerca de 1 minuto a correr 
tic

%Constantes 
N=5e4; 
xc=-30;
dt=10;
tf=50000;
D=1/4; 

t=0:tf;
F_c=zeros(1,length(t));
S_c=F_c; 
for n=1:N
    r=zeros(length(t),2); 
    x=0;
    y=0;
    survival=1;
    for i=2:length(t)  %O mesmo ciclo for do ex anterior  
        A=randi([1,4]);  
        mov_x=0;
        mov_y=0;
        if A==1
            mov_y=1;
        elseif A==2
            mov_x=1;
        elseif A==3
            mov_y=-1;
        elseif A==4
            mov_x=-1;
        end
        x=mov_x+x; 
        y=mov_y+y;
        %Detetar traps
        if x==xc
            F_c(i)=F_c(i)+1;
            break
        else
            S_c(i)=S_c(i)+1;
        end
    end    
end

t_dt=1:dt:tf-2;
idx=1;
for i= 0:dt:length(t)-dt
    s1=0;
    s2=0;
    for j=1:dt
        s1 = s1 + F_c(i+j);
        s2 = s2 + S_c(i+j);
    end
    N_f(idx)=s1;
    N_s(idx)=s2;
    idx=idx+1;
end


%Normalizar
F=N_f./(N*dt);
figure()
subplot(1,2,1)
plot(t_dt,F,'.')
hold on

%Valores teóricos de F
F_anal=abs(xc)./(sqrt(4*pi*D*t_dt.^3)) .* exp(-((xc)^2)./(4*D.*t_dt));

plot(t_dt,F_anal,'k-','LineWidth',1.5)
title('Probabilidade da Primeira Passagem, F(t)')
legend('F','F_{teórico}','Location','best')
xlabel('t'),ylabel('F')
subplot(1,2,2)
plot(log10(t_dt),log10(F),'.')
hold on
plot(log10(t_dt),log10(F_anal),'k-','LineWidth',1.5)
title('Probabilidade da Primeira Passagem, F(t) (log-log)')
legend('log_{10}(F)','log_{10}(F_{teórico})')
xlabel('log_{10}t'),ylabel('log_{10}F')

%Normalizar
S=(N_s./(N*dt));

figure

subplot(1,2,1)
plot(t_dt(2:end),S(2:end),'b--','LineWidth',3)
hold on
S_anal= erf(30./(2*sqrt(D.*t_dt)));
plot(t_dt,S_anal,'k-','LineWidth',1.5)
title('Probabilidade de Sobrevivência, S(t)')
legend('S','S_{teórico}','Location','best')

subplot(1,2,2)
plot(t_dt(2:end),S(2:end),'r--','LineWidth',3)
hold on
%Valores teóricos de S
S_anal= erf(30./(2*sqrt(D.*t_dt)));
plot(t_dt,S_anal,'k-','LineWidth',1.5)
title('Probabilidade de Sobrevivência, S(t) [com zoom]')
legend('S','S_{teórico}','Location','best')
xlim([2500 2700])

toc
