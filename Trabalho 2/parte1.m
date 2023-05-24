close all
clear all
clc

%constantes
p=0.5;
q=p;
ti=0;
tf=40;


rngs=[rng(0),rng(1),rng(2)];

N=5000000;

%Definir Vetores
t=ti:tf;
x(1)=0;



%Definir x(t)

weigths=[p q];
population=[1 -1];


for j=1:length(rngs)
    rng(j)
    for i=2:length(t)
        S=randsample(population,N,true,weigths);
        x(i,j)=sum(S(1:i));
    end
end
plot(t,x,'-o')
title('1D Random Walks')
xlabel('Tempo')
ylabel('Posição em x')
legend(['Trajecto 1'],['Trajecto 2'],['Trajecto 3'],'Location','Southwest')




%% Task 1.2
clc
clear
close all

N=50000; 
x=-500:500;


matrz_tf=[40,41;400,401;4000,4001];
p = 0.5; 


P_media=zeros(length(x),3);
for i=1:3
    P=zeros(length(x),2);
    lista_tf=matrz_tf(i,:); 
    for n=1:N
        for j=1:length(lista_tf)
        tf=lista_tf(j);
        t=0:tf;
        posicao=zeros(length(t),1);
            for k=1:length(t)-1   
                R=rand;   
                if R<p
                    S=-1;
                else
                    S=1;
                end
            posicao(k+1)=S+posicao(k);
                if k==length(t)-1
                    xi=posicao(k+1);
                    idx=xi+abs(min(x))+1;
                    P(idx,j)=P(idx,j)+1;
                end
            end
        end
    end

    for k=1:length(P)
        P_media(k,i)=(P(k,1)+P(k,2))/2;
    end
    P_media(:,i)=P_media(:,i)./N;
end
figure(1)
plot(x,P_media(:,1),'b',x,P_media(:,2),'r',x,P_media(:,3),'y','LineWidth',1.7)
title('<P(x,t)>')
legend('t = 40,41' ,'t = 400,401' ,'t = 4000,4001')
xlabel('x')
ylabel('Probabilidades')
xlim([-200 200])

tf=(lista_tf(1)+lista_tf(2))/2;
P_t=1/sqrt(2*pi*tf).*exp(-(x.^2)/(2.*tf));

figure(2)
plot(x,P_media(:,3),x,P_t,'-','LineWidth',1.7)
title('Valores Esperados de P(x,t) and <P(x,t)>')
legend('<P(x,t)>','P_{t}')
xlabel('x')
ylabel('Probabilidades')
xlim([-200,200])