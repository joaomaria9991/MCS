close all
clear variables
clc

H = 0.1;    %Campo magnético
J = 1;      
N = 1000;   
m = 100000; 

Temp = 0.1:0.05:10;

M = zeros(m,1);
MT = zeros(length(Temp),1);

S = ones(N,1);



for k = 1:length(Temp)
   M100 = zeros(length(Temp),m/100-200);
  for i = 1:m
        p_centro = randi(N);

          %Garantir o anel
        if p_centro > 1 && p_centro < N
            p_previous = p_centro - 1;
            p_next = p_centro + 1;
        elseif p_centro == 1
            p_previous = N;
            p_next = p_centro + 1;
        else
            p_previous = p_centro - 1;
            p_next = 1;
        end

        Spin        = S(p_centro);
        SpinLeft    = S(p_previous);
        SpinRight   = S(p_next);

        %%Virar S
        E      = -J*Spin*(SpinLeft + SpinRight) - H*Spin;
        E_new   = -J*(-Spin)*(SpinLeft + SpinRight) - H*(-Spin);
        dE = E_new - E;

        if dE <= 0
            S(p_centro) = -S(p_centro);
        else  
            beta = 1/Temp(k);
            w = exp(-beta*dE);
            r = rand(1);
            if r <= w
                S(p_centro) = -S(p_centro);
            end
        end

        M(i) = (1/N)*sum(S);
  end
    % Magnetização média
    l=1;
    for i = 100:100:m-100  %Escolher indices pretendidos
        M100(k,l) = (1/i)*sum(M(1:i+99));
        l = l + 1;
    end
    
    MT(k) = mean(M100(k,:));
end


M_theo = zeros(length(Temp),1);
for j = 1:length(Temp)
    beta = 1/Temp(j);
    M_theo(j) = (sinh(beta*H)) / ( sqrt( (sinh(beta*H))^2 + exp(-4*beta*J)) );
end

%% Plotting
close all


figure
plot(Temp,M_theo,'k-');
hold on, grid on
plot(Temp,MT,'r.-');
xlabel('Temperatura'); ylabel('Magnetização');
title('Magnetização em função Temperatura')
legend('M_{Teórica}','M_{Computada}')






