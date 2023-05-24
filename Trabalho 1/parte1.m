close all
clear all
clc


N=[10^2 10^4 10^6 ];


for j=1:length(N)

M=100;
N_val=N(j);

for i=1:N_val
    x(i)=randi(M);
end


%media
med=sum(x)/N_val;

med_anal=(M+1)/(2);

med_conv(j)=abs(med-med_anal);



%variancia

for i=1:N_val
    dif(i)=(x(i)-med)^2;
end

sigma=(1/N_val)*sum(dif);
sigma_anal=(M^2-1)/(12);
sigma_conv(j)=abs(sigma-sigma_anal);


%Determinar p

y=x(x<60);
p=length(y)/N_val;
p_anal=0.6;
p_con(j)=abs(p-p_anal);

end



