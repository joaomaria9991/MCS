close all
clear all
clc



N_val=[10^2,10^4,10^6];


for i=1:length(N_val)
    
    N=N_val(i);

    x=rand(1,N);
    y=rand(1,N);
    z=x.*y;

    med_x=(sum(x)/N);
    med_y=(sum(y)/N);
    med_z=(sum(z)/N);


    med_z_anal=med_x.*med_y;

    conv(i)=abs(med_z_anal-med_z);

end