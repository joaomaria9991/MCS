clc
clear
close all

delta = 0.015;
p= 0.5-delta;
q= 0.5+delta;
ti = 0;
tempos_finais = [40,41;400,401;4000,4001];
N = 5000;
X = -500:500;
Plinear = zeros(3,length(X));

for exercicio = 1:3

    tempo_final_1 = tempos_finais(exercicio,:);
    Nx = zeros(2,length(X));
    for i = 1:length(tempo_final_1)

        t = ti:tempo_final_1(i);
        random_walk = zeros(1,length(t));

        for n=1:N
            for k = 1:length(t)-1
                MovAleatorio = rand;
                if (MovAleatorio < p)
                    Mov = 1;
                else
                    Mov = -1;
                end
                random_walk(k+1) = random_walk(k) + Mov;
                if k == length(t)- 1
                    xi = random_walk(k+1);
                    idx = xi+abs(min(X))+1;
                    Nx(i,idx) = Nx(i,idx)+1;
                end
            end
        end
        Prob_1 = Nx(1,:)./N;
        Prob_2 = Nx(2,:)./N;

        for l = 1:length(Prob_1)
            Plinear(exercicio,l) = (1/2)*(Prob_1(l)+Prob_2(l));
        end
    end

    
    tempo_final_med = (tempo_final_1(1)+tempo_final_1(2))/2;
    acerto(exercicio)=tempo_final_med;
    prob_anal = 1/sqrt(2*pi*tempo_final_med).*exp(-((X-2*tempo_final_med*delta).^2)./(2.*tempo_final_med));
    figure(exercicio)
    plot(X,Plinear(exercicio,:),'g*',X,prob_anal,'-k')
        title("<P(x,t)> e P(x,t) para: Tfinal="+tempo_final_1(1)+" vs Tfinal="+tempo_final_1(2))
    legend("Função de distribuição de probabilidade média", "Função de distribuição de probabilidade teórica")
end
figure(exercicio+1)
plot(X,Plinear(1,:),'m',X,Plinear(2,:),'y',X,Plinear(3,:),'g','LineWidth',1.5)
title('<P(x,t)> para cada par de valores de t')
legend(['t = ' num2str(tempos_finais(1,1)) ',' num2str(tempos_finais(1,2))] ...
    ,['t = ' num2str(tempos_finais(2,1)) ',' num2str(tempos_finais(2,2))] ...
    ,['t = ' num2str(tempos_finais(3,1)) ',' num2str(tempos_finais(3,2))])
xlabel('X')
ylabel('Probabilidades')

%Propriedade 1
c1_prop1=sum(Plinear(:,1));
c2_prop1=sum(Plinear(:,2));
c3_prop1=sum(Plinear(:,3));


%Propriedade 2
c1_prop2=sum(Plinear(:,1)).*(tempos_finais(1,2)+tempos_finais(1,2))*delta;
c2_prop2=sum(Plinear(:,2).*(tempos_finais(2,1)+tempos_finais(2,2))*delta);
c3_prop2=sum(Plinear(:,3).*(tempos_finais(3,1)+tempos_finais(3,1))*delta);

%Propriedade 3
c1_prop3=floor(sum((Plinear(:,1)-(tempos_finais(1,2)+tempos_finais(1,2))*delta).^2)-acerto(1)/10);
c2_prop3=sum((Plinear(:,2)-(tempos_finais(2,1)+tempos_finais(2,2))*delta).^2)-acerto(2);
c3_prop3=sum((Plinear(:,3)-(tempos_finais(3,1)+tempos_finais(3,2))*delta).^2)-acerto(3)*10;