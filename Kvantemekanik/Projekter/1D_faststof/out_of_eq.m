%her beregner vi osc. som følge af at vi bringer atomnummer l ud af
%ligevægt med længde b.
clear all; close all; clc;


%antal atomer: 
N = 50; 
l = 1;

%afstand mellem atomer ved hvilelængde:
a = 1;

%afstand fra ligevægt:
b = 2/3*a; 

%k-værdier:
k = [1:N/2]*2*pi/(N*a);

k_even = [2:2:N/2]*2*pi/(N*a);

%grundfrekvens:
W = 2;

%dispersionsrelation:
w = 2*W*sin(k*a/2); 

w_even = 2*W*sin(k_even*a/2);


%tider:
t = linspace(0, 60/W, 1e4)';
 
%farter:
v_phase = w./k;
v_small_k = W*a;
v_group = W*a*cos(k*a/2);


x = zeros(N,length(t));
inh = 2; 
for j = 1:N
    for m = 1:length(t)
        x(j,m) = inh*2*b/N*sum(cos(k*(j-l)*a).*cos(w*t(m)))+j*a;
    end
end

% x_half_uneq = zeros(N,length(t));
% inh_2 = 2; 
% for j = 1:N
%     for m = 1:length(t)
%         x_half_uneq(j,m) = -inh_2*2*b/N*sum(cos(k_even*j*a-w_even*t(m)))+j*a;
%     end
% end


% figure 
% hold on
% set(gca,'fontsize', 12)
% xlabel('\Omega t')
% ylabel([num2str(inh) '<x>'])
% 
% for j = 1:N
%     plot(W*t,v_small_k*t + l*a - ((j-1)*N-1)*a ,'r-','linewidth',0.8)
%     %plot(W*t,(-v_small_k*t + (j*N+1)*a),'b-','linewidth',0.8)
%     plot(W*t,x(j,:),'k-','linewidth',1.2)
%     %legend('Forward wave','Backward wave','Ground state osc.')
% end
% axis([0 W*max(t) 0 (N+2)*a])


figure 
hold on
set(gca,'fontsize', 12)
xlabel('\Omega t')
ylabel([num2str(inh) '<x>'])

for j = 1:N
    plot(W*t,v_small_k*t  - ((j-1)*N-1)*a ,'r-','linewidth',0.8)
    plot(W*t,(-v_small_k*t + (j*N+1)*a),'b-','linewidth',0.8)
    plot(W*t,x(j,:),'k-','linewidth',1.2)
    legend('Forward wave','Backward wave','Ground state osc.')
end
axis([0 W*max(t) 0 (N+2)*a])