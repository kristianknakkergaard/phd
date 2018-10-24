clear all; close all; clc;

%antal atomer: 
N = 50; 

%afstand mellem atomer ved hvilelængde:
a = 1;


%oscillationer relativt til atomnummer:
l = 1;

%k-værdier:
k = [1:N/2]*2*pi/(N*a);

%grundfrekvens:
W = 2;

%dispersionsrelation:
w = 2*W*sin(k*a/2); 


%tider:
t = linspace(0, 2e2/W, 1e5)';

x_gs1 = zeros(length(k),length(t));
x_gs2 = zeros(length(k),length(t));

n_prime = 10;
x_exc = x_gs1; 

gs_inh = 5;
exc_inh = 5;

for j = 1:N
    x_gs1(j,:) = gs_inh*1/N*sum(cos(kron(w,t)).*kron(1./w.*cos(k*(j-l)*a),ones(length(t),1)) , 2) + j*a;
    x_gs2(j,:) = gs_inh*1/N*sum(sin(kron(w,t)).*kron(1./w.*cos(k*(j-l)*a),ones(length(t),1)) , 2) + j*a;
    x_exc(j,:)= exc_inh*1/(N*w(n_prime))*cos(k(n_prime)*(j-l)*a-w(n_prime)*t);
end

x_exc_2 = x_exc + x_gs1;

v_phase = w./k;
v_small_k = W*a;
v_group = W*a*cos(k*a/2);
% 
figure 
hold on
set(gca,'fontsize', 12)
xlabel('\Omega t')
ylabel([num2str(gs_inh) 'Cov_0(x_j,x_1)'])

for j = 1:N
    plot(W*t,(v_small_k*t - ((j-1)*N-1)*a),'r-','linewidth',0.8)
    plot(W*t,(-v_small_k*t + (j*N+1)*a),'b-','linewidth',0.8)
    plot(W*t,x_gs1(j,:),'k-','linewidth',1.2)
    legend('Forward wave','Backward wave','Ground state osc.')
end
axis([0 W*max(t) 0 (N+1)*a])

figure 
hold on
set(gca,'fontsize', 12)
xlabel('\Omega t')
ylabel([num2str(gs_inh) 'Cov_0(x_j,x_1)'])

for j = 1:N
    plot(W*t,(v_small_k*t - ((j-1)*N-1)*a),'r-','linewidth',0.8)
    plot(W*t,(-v_small_k*t + (j*N+1)*a),'b-','linewidth',0.8)
    plot(W*t,x_gs2(j,:),'k-','linewidth',1.2)
    legend('Forward wave','Backward wave','Ground state osc.')
end
axis([0 W*max(t) 0 (N+1)*a])


% figure 
% hold on
% set(gca,'fontsize', 12)
% xlabel('\Omega t')
% ylabel('5 x_{exc}')
% plot(W*t,w(n_prime)/k(n_prime)*t+a,'r-')
% for j = 1:N
%     plot(W*t,x_exc(j,:)+j*a,'k-')
% end
% legend('Traveling wave','Excited osc.')
% axis([0 W*max(t) 0 (N+2)*a])
% 
% figure 
% hold on
% for j = 1:N
%     plot(t,x_exc_2(j,:),'b-')
% end

%her beregner vi osc. som følge af at vi bringer atomnummer l ud af
%ligevægt med længde b.

% l = N; 
% b = 2/3*a; 
% 
% x = zeros(N,length(t));
% inh = 3; 
% for j = 1:N
%     for m = 1:length(t)
%         x(j,m) = inh*4*b/N*sum(cos(k*j*a-w*t(m)).*cos(k*l*a))+j*a;
%     end
% end
% 
% 
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

% m_prime = 16;
% j_prime = 35; 
% 
% filter = t<20; 
% figure
% hold on
% plot(W*t(filter),x(j_prime,filter), 'k-')
% plot(W*t(filter),0.1*cos(k(m_prime)*(j_prime*a-v_small_k*t(filter)) - 0.35 )+j_prime*a-0.085, 'r-')
% 
% 


