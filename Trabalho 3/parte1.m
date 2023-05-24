%Projeto 3

%% Task 1.1
close all
clear all
clc


tf=100;
n=3;

t=0:tf;
r=zeros(length(t),2,n);


coluna=['g','b','r','y'];


for k=1:n
    x=0;y=0;
     for i=2:length(t)
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
            r(i,:,k)=[x y]; % posição da partícula
     end


    figure(1)
    plot3(t,r(:,1,k),r(:,2,k),[coluna(k) '-o'])
    hold on
    figure(2)
    plot(t,r(:,1,k),[coluna(k) '-o'])
    hold on
    figure(3)
    plot(t,r(:,2,k),[coluna(k) '-o'])
    hold on
    figure(4)
    plot(r(:,1,k),r(:,2,k),[coluna(k) '-o'])
    hold on
end

figure(1)
title('Posição em função do tempo')
legend('Passeio Aleatório 1','Passeio Aleatório 2','Passeio Aleatório 3','Location','best')
xlabel('t'),ylabel('x'),zlabel('y')
figure(2)
title('Posição em função do tempo (No plano Y=0)')
legend('Passeio Aleatório 1','Passeio Aleatório 2','Passeio Aleatório 3','Location','best')
xlabel('t'),ylabel('x')
figure(3)
title('Posição em função do tempo (No plano X=0)')
legend('Passeio Aleatório 1','Passeio Aleatório 2','Passeio Aleatório 3','Location','best')
xlabel('t'),ylabel('y')

figure(4)
title('Trajetorias dos passeios')
legend('Passeio Aleatório 1','Passeio Aleatório 2','Passeio Aleatório 3','Location','best')
xlabel('x'),ylabel('y')




%% Task 1.2
close all
clear all
clc


tf_list=[5000, 5001];
xx=-500:500;
yy=xx;

N=10000; %nr of repetitions

Nxy=zeros(length(xx),length(yy),2);
Pmed=zeros(length(xx),length(yy));
s1=0;
s2=0;
s3=0;
s4=0;


for n=1:N
    for j=1:length(tf_list)
        tf=tf_list(j);
        t=0:tf;
        x=0;
        y=0;
        for i=2:length(t)
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
        end
        id_x=x+abs(min(xx));
        id_y=y+abs(min(yy));
        Nxy(id_x,id_y,j)=Nxy(id_x,id_y,j)+1;
    end
end

P_anal= zeros(length(xx),length(yy));
t_anal=mean(tf_list);
c=1/(pi*t_anal);



for i=1:length(xx)
    for j=1:length(yy)
        P_anal(i,j)=c*exp(-((xx(i)^2+yy(j)^2))/t_anal);
        Nmean=(Nxy(i,j,1)+Nxy(i,j,2))./2;
        if Nmean==0
            Pmed(i,j)=NaN;
        else
            Pmed(i,j)=Nmean./N;
            s1=s1+Pmed(i,j);
            s2=s2+Pmed(i,j)*xx(i);
            s3=s3+Pmed(i,j)*yy(i);
            s4=s4+Pmed(i,j)*(xx(i)^2+yy(i)^2);
        end
    end
end



[x,y]=meshgrid(xx,yy);
figure
subplot(2,2,1)
plot3(x,y,Pmed,'ko')
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(2,2,2)
mesh(x,y,P_anal)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(2,2,3)
plot3(x,y,Pmed,'ko'),view(0,0)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(2,2,4)
mesh(x,y,P_anal),view(0,0)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
sgtitle('Distribuição de probabilidade no ponto (x,y)')

%verificar as propriedades pedidas
s1
s2
s3
s4

