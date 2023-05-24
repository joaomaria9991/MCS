
%%  


close all
clear 
clc

%constantes
N = 10000;
c = 50;
m = 100;
L = (N/2)*c;


A = zeros(N,N);
q = zeros(1,N);


for k = 1:L
    i = 0;
    j = 0;
    while i == j || A(i,j) == 1
        i = randi([1,N]);
        j = randi([1,N]);
    end
    
    A(j,i) = 1;
    A(i,j) = 1;
end

for l = 1:N
    q(l) = sum(A(l,:));
end

q_med = (1/N)*sum(q)
B = (1/(N*q_med))*sum(q.*(q-1))
q_med2 = (1/N)*sum(q.^2)
q_med3 = (1/N)*sum(q.^3)

it = 1:100;


for qi = it
    Pq(qi) = P(qi,q);
    Pq_anal(qi) = exp(-c)*((c^qi)/factorial(qi));
end

q_teo = c;
B_teo = q_teo;
const1 = B/q_med


plot(it,Pq,'gx');
hold on
plot(it,Pq_anal,'-k');
xlabel('q')
ylabel('P')
title('Distribuição de probabilidade de graus')
legend('P(q)_{Experimental}','P(q)_{teórica}','Location','best')

%%  
npt = (1/6)*N*((1/N)*sum(q.*(q-1)));

ntr = (1/6)*(q_med^3); %melhor que trace

C = ntr/npt;
C_teo = q_med/N;
const2 = C/C_teo


%% 
Q = q_med2/q_med;
sigma_sq = q_med3/q_med - (q_med2/q_med)^2;

rho = 0;
for i = 1:N
    for j = 1:N
        rho = rho + (A(i,j)*(q(i)-Q)*(q(j)-Q));
    end
end
rho = rho/(N*q_med*sigma_sq)  % normalização -1 a 1




function [Pq] = P(g, vetor) % grau do vertice
    Nq = sum(vetor==g);
    N = length(vetor);
    Pq = Nq/N;
end