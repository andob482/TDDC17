close all; clear all; 
%%
%%Low-pass low degree
N = 2^10;
x = randn(1, N);                               %vitt brus 
X = fft(x, N);                                 %bruset i frekvensdom�n
theta_axis = linspace(0,1,N);                  %for plotting
R0 = 1;
time_axis = linspace(0,1,N);
H_z_low = 0.5*(1+exp(-1i*2*pi*theta_axis));    %low degree low_pass theoretical
Y_z_low = X.*H_z_low;                          %utsignal frekvensdom�n
y_t_low = ifft(Y_z_low);

Ry_theo_low = 0.5*(1 + cos(2*pi.*theta_axis)); %PSD theoretical low deg
ry_theo_low = ifft(Ry_theo_low);                   %ACF theoretical low deg


Ry_est = abs(ACF_est(y_t_low));
ry_est = ifft


figure(1)
plot(theta_axis, Ry_theo_low)                  %Plotting Ry_low deg
figure(2)                                    
plot(time_axis, fftshift(ry_theo_low))
figure(3)
plot(theta_axis, fftshift(Ry_est))

%Plotting ry_low deg
%%
%Lowpass high degree

[b, a] = butter(25, 0.1);

y_hd_lp = filter(b, a, x);

figure(2)
plot(time_axis, y_hd_lp)



