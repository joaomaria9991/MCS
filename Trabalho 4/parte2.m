clear
clc

%Constantes
N = 100;      
num_mstate = 5000;     
dT = 0.05;      
Temp = 0.1:dT:10;  
num_T = length(Temp);
H = [0, 0.001];       
num_field = length(H);

%Inicializar vetores
m_comp = zeros(num_T, num_field);
chi_exp = zeros(num_T, num_field);
m_anal = zeros(num_T, num_field);
chi_anal = zeros(num_T, num_field);

idx_Plot = fix(num_T / 5) * (1:5);

J = 1;     

for h = 1:num_field
    for i = 1:length(Temp)
        m_avg = zeros(num_mstate / 100, 1);
        x_avg = zeros(num_mstate / 100, 1);
        m = zeros(num_mstate, 1);
        sigma = ones(N, 1);      % Generate all spins up

        for j = 1:num_mstate
            estado = randi(N);   % Choose random spin
    
            E = getEnergyLR(sigma, estado, H(h), J, N);
    
            sigma(estado) = -sigma(estado);     % Flip spin
            
            E_new = getEnergyLR(sigma, estado, H(h), J, N);
    
            %Algoritmo de Metropolis
            dE = E_new - E;
    
            if dE > 0
                w = exp(-dE / Temp(i));
                r = rand();
                if r >= w
                    sigma(estado) = -sigma(estado);
                end
            end
            m(j) = (1/N) * sum(sigma);
            if mod(j, 100) == 0 && j > 100
                m_avg(j/100 - 1) = (1 / j) * sum(m(1:j));
                x_avg(j/100 - 1) = (N / Temp(i)) * ((1 / j) * sum(m(1:j) .^ 2) - m_avg(j/100 - 1) ^ 2);
            end
        end  
        m_comp(i, h) = m_avg(end - 1);
        chi_exp(i, h) = x_avg(end - 1);    
        % Valores teoricos, analitico
        syms m;
        ftmp = (m - tanh(((J * m) / Temp(i)) + (H(h) / Temp(i)))) == 0;
    
        m_anal(i, h) = vpasolve(ftmp, m, 1);
        chi_anal(i, h) = 1 / (Temp(i) * cosh(((J * m_anal(i, h)) / Temp(i)) + (H(h) / Temp(i))) ^ 2 - (J / Temp(i)));
    end
end


%% Plots

close all

figure(1)
plot(Temp, chi_exp(:, 1), 'r.-', 'DisplayName', 'Valores Computados', 'LineWidth', 1.2),  hold on
plot(Temp, chi_anal(:, 1), 'k-', 'DisplayName', 'Valores Teóricos', 'LineWidth', 1.2)
xline(J, 'g--', 'DisplayName', 'Temperatura Crítica = J', 'LineWidth', 1.5)
xlabel('Temperatura'), ylabel('Suscetibilidade Média')
title('Suscetibilidade em função Temperatura (H = 0)')
legend()
xticks(0:10)

figure(2)
plot(Temp, m_comp(:, 2), 'r.-', 'DisplayName', 'Valores Computados', 'LineWidth', 1.2), hold on
plot(Temp, m_anal(:, 2), 'k-', 'DisplayName', 'Valores Teóricos', 'LineWidth', 1.2)
xline(J, 'g--', 'DisplayName', 'T_c = J', 'LineWidth', 1.5)
xlabel('Temperatura'), ylabel('Magnetização')
title('Magnetização em função Temperatura (H = 0.001)')
xticks(0:10)

figure(3)
plot(Temp, m_comp(:, 2), 'r.-', 'DisplayName', 'Valores Computados', 'LineWidth', 1.2), hold on
plot(Temp, m_anal(:, 2), 'k-', 'DisplayName', 'Valores Teóricos', 'LineWidth', 1.2)
xline(J, 'g--', 'DisplayName', 'T_c = J', 'LineWidth', 1.5)
xlabel('Temperatura'), ylabel('Magnetização')
xlim([0.5 1.5])
title('Magnetização em função Temperatura (H = 0.001)')

figure(4)
plot(Temp, chi_exp(:, 2), 'r.-', 'DisplayName', 'Experimental', 'LineWidth', 1.2), hold on
plot(Temp, chi_anal(:, 2), 'k-', 'DisplayName', 'Theoretical', 'LineWidth', 1.2)
xline(J, 'g--', 'DisplayName', 'T_c = J', 'LineWidth', 1.5)
xlabel('Temperatura'), ylabel('Suscetibilidade Média')
xlim([0.5 1.5])
title('Suscetibilidade em função Temperatura (H = 0.001)')



%% Funcao
function E = getEnergyLR(sigma, idx, H, J, n)
    sigma_sum = sum(sigma) - sigma(idx);
    E = -(J/n) * sigma(idx) * sigma_sum - H * sigma(idx);
end