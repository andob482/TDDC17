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
%%%%%%%%%%%Lowpass high degree%%%%%%%%%%%%%%
N = 2^10;
%Insignal
wg_noise = randn(1, N);                               %vitt brus 
wg_noisefreq = fft(wg_noise, N);                      %bruset i frekvensdomän

%Till plott
theta_axis = linspace(0,1,N);                         %for plotting
R0 = 1;
time_axis = linspace(0,1,N);
Ts = 1;                                                 %Step for the vectos
nn = ((-N)/2)+1:Ts:(N)/2;                               %Integer vector for plotting
ff= linspace(0,1,N);
k_axis = (-N+1:N-1); %Måste ha detta pga dubbelsidig ACF
freq_axis_dubble = linspace(0,1,2*N-1); %Måste ha detta pga dubbelsidig ACF

%High degree filter 
[b, a] = butter(25, 0.4);  %Butterfilter med cut off 0.2
                

%%%%%%%%%Filtered output
%Tidsdomän
y_hd_lp = filter(b, a, wg_noise);           %filtrerad utsignal tidsdomän
%Frekvensdomän
Y_hd_lp = fft(y_hd_lp);                        %filrerad Hz

%%%%%%%%%ESTIMERING
%ACF
ACF_est = ACF_two_sided(y_hd_lp, N);

%PSD
PSD_est = abs(fft(ACF_est));
freq_axis_dubble = linspace(0,1,2*N-1);
%ACF
%ry_hd_lp= ifft(Ry_hd_lp);

%ry_hd_lp_est = ACF_est(y_hd_lp, N);
%Ry_hd_lp_est = abs(fft(ry_hd_lp_est));

%%%%%%%%%TEORETISK
%ACF_ter
theta_c = 2*0.2; %ger cutoff på 0,2
ACF_ter = theta_c*sinc(theta_c*nn);

%PSD_ter
PSD_ter = abs(fft(ACF_ter));

%figure(2)
%plot(nn, abs(y_hd_lp))
%figure(3)
%plot(ff, Ry_hd_lp)

%%%%%%%%%%%Plottar HIGH DEGRE
%ACF_est PLOT
figure(4)
stem(k_axis, ACF_est)
xlim([-20, 20])
title('Estimerad ACF');

%PSD_est PLOT
figure(5)
plot(freq_axis_dubble, PSD_est)
title('Estimerad PSD')

%ACF_ter PLOT
figure(6)
stem(nn, ACF_ter)
xlim([-20, 20])

figure(7) 
plot(ff, PSD_ter)


