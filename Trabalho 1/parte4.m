close all
clear all
clc

%Parametros forncidos
M=9;
N=21;
k=[10^2 10^4 10^6];
n_total=0:N;


for i=1:length(k)
    for j=1:k(i)
        x=randi(M,[N,1]);
        n(j)=numel(x(x==3));

    end
    p=1/M;
    for h=0:N
        Ntr=sum((n==h)==1); %dupla indexação lógica
        P(h+1)=Ntr/k(i);
        Bin(h+1)=nchoosek(N,h)*(p^h)*(1-p)^(N-h); %Binomial
        Poi(h+1)=exp(-N*p)*((N*p).^h)/factorial(h); %Poisson
        Gau(h+1)=(1/sqrt(2*pi*N*p))*exp(-((h-N*p)^2)/(2*N*p)); %Gaussiana


    end

        figure(i)
        plot(n_total,P,'-b')
        hold on
        plot(Bin,'-r')
        plot(Poi,'-g')
        plot(Gau,'-k')
        legend("Probabilidade","Binomial","Poisson","Gaussiana")



end