close all
clear all
clc



N=[10^2,10^4,10^6];

dy=0.005;
y=0:dy:1;


for N_vals=1:length(N)

    N_sim=N(N_vals);

   x=sqrt(rand(N_sim,1));
   %x=rand(N_sim,1);



    for k=1:length(y)-1

        lower=y(k);
        upper=y(k+1);
        p(k,N_vals)=length(x(x<=upper & x>=lower));
        bin_center(k)=(upper+lower)/2;

        

    end



    p(:,N_vals)=p(:,N_vals)/(sum(p(:,N_vals))*(dy));

    m(:,N_vals)=sum(bin_center*p(:,N_vals)*dy);

    sig(:,N_vals)=sum((bin_center(:)-m(:,N_vals)).^(2).*p(:,N_vals))*dy;

end



subplot(1,3,1)
plot(bin_center,p(:,1))
title('N=100')
hold on

subplot(1,3,2)
plot(bin_center,p(:,2))
title('N=10^4')
ylim([0,3])
hold on

subplot(1,3,3)
plot(bin_center,p(:,3))
title('N=10^4')
ylim([0,3])











