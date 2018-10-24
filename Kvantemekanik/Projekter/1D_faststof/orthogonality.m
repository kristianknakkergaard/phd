close all; clc;
N = 6; 
n = 3; 
n_prime = 1; 

j = [1:N];

f_real = cos((n+n_prime)*2*pi/N*j);
f_imag = sin((n+n_prime)*2*pi/N*j);


f_real_sum = sum(f_real)
f_imag_sum = sum(f_imag)

figure
plot(j, f_real,'r.')

figure 
plot(j, f_imag,'k.')