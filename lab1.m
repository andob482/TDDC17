close all; clear all; 

%%Low-pass low degree
N = 2^10;
wg_noise = randn(1, N);                               %vitt brus 
wg_noisefreq = fft(wg_noise, N);                      %bruset i frekvensdomï¿½n
theta_axis = linspace(0,1,N);                         %for plotting
R0 = 1;
time_axis = linspace(0,1,N);

Ts = 1;                               %Step for the vectos
nn = ((-N)/2)+1:Ts:(N)/2;             %Integer vector for plotting
ff= linspace(0,1,N);                  %Frequency vector for plotting

H_z_low = 0.5*(1+exp(-1i*2*pi*theta_axis));           %low degree low_pass theoretical
Y_z_low = wg_noisefreq.*H_z_low;                      %utsignal frekvensdomï¿½n
y_t_low = ifft(Y_z_low);

Ry_theo_low = 0.5*(1 + cos(2*pi.*theta_axis));      %PSD theoretical low deg
ry_theo_low = ifft(Ry_theo_low);                    %ACF theoretical low deg


ry_est = ACF_est(y_t_low, N);
Ry_est = abs(fft(ry_est));


figure(1)                            %Plotting Ry_low deg
subplot(2,1,1)
plot(ff, Ry_theo_low)
title('Teoretiskt PSD')
xlabel('\theta')
ylabel('Ry(\theta)')          
subplot(2,1,2)
plot(ff, Ry_est)
title('Estimerad PSD')
xlabel('\theta')
ylabel('Ry(\theta)') 

figure(2)                           %Plottar av teoretiskt och estimerat ACF
subplot(3,1,1)                                    
stem(nn, fftshift(ry_theo_low)) 
xlim([-20,20])
title('Teoretiskt ACF')
xlabel('K')
ylabel('rx[k]')
subplot(3,1,2)
plot(nn, ry_est)
title('Estimerad ACF')
xlabel('k')
ylabel('ry[k]') 
subplot(3,1,3)
stem(nn,abs(ry_est))
xlim([-20,20])
title('Estimerad ACF')
xlabel('k')
ylabel('ry[k]') 

%%
clear all; close all;
%Lowpass high degree
N = 2^10;
wg_noise = randn(1, N);                               %vitt brus 
wg_noisefreq = fft(wg_noise, N);                      %bruset i frekvensdomï¿½n
theta_axis = linspace(0,1,N);                         %for plotting
R0 = 1;
time_axis = linspace(0,1,N);

Ts = 1;                                                 %Step for the vectos
nn = ((-N)/2)+1:Ts:(N)/2;                               %Integer vector for plotting
ff= linspace(0,1,N); 

[b, a] = butter(25, 0.2);                       %filter
H_hd_lp = butter(25,0.2);                       %filter?

Y_hd_lp = filter(b, a, wg_noisefreq);           %filtrerad utsignal frekvensdomän
y_hd_lp = ifft(Y_hd_lp);                        %filrerad utsignal tidsdomän

Ry_hd_lp = abs(H_hd_lp).^2;
ry_hd_lp= ifft(Ry_hd_lp);

ry_hd_lp_est = ACF_est(y_hd_lp, N);
Ry_hd_lp_est = abs(fft(ry_hd_lp_est));

%figure(2)
%plot(nn, abs(y_hd_lp))
figure(3)
plot(ff, Ry_hd_lp)
figure(4)
plot(nn, ry_hd_lp_est)
figure(5)
plot(ff, Ry_hd_lp_est)



