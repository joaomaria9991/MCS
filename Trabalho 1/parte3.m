close all
clear all
clc


n=[10,100, 1000];
N=1e6;
dy=0.005;
k=dy:dy:1;
y_centre_bins=k+(0.5*dy);


Y_med=zeros(1,length(n));
x_med=zeros(1,length(n));
x_sig=zeros(1,length(n));


for i=1:length(n)


    for j=1:N
        x=rand(1,n(i));
        Y_med(j)=(1/n(i))*sum(x);
    end

    [M,vert]=histcounts(Y_med,k);
    P=M/(N*dy);

    norm=sum(P*dy);

    figure(i)
    plot(y_centre_bins(1:end-1),P)
    

    y_mean(i)=sum(vert(1:end-1).*P*dy);
    y_sig(i)=sum(dy*P.*(vert(1:end-1)-y_mean(i)).^2);


    x_med(i)=(1/N)*sum((1/n(i))*sum(x));
    x_sig(i)=(1/N)*sum((1/n(i))*sum((x-x_med(i)).^2));


    conv(i)=abs(y_sig(i)-x_sig(i)); %Converge para zero


end


