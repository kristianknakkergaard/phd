% ------------------------------------------------------------------
%  Denne Matlab-kode skal illustrere forskellige løsninger til den delta-
%  funktions potentialet som beskrevet i kapitel 2.5 i Griffith's
%  "Introduction to Quantum Mechanics".
%
%  Bemærk, vi har sat hbar = m = 1 i de aktuelle formler.
% ------------------------------------------------------------------
close all
clear all
FontSize = 14;

% --- Styrke af delta-funktions-potentialet [V(x) = -alpha*delta(x)]:
alpha = 2;

% --- Definer først en x-akse: ---
xmin = -20;
xmax = 20;
dx = 0.02;
x = [xmin:dx:xmax];

% --- Definer så en start-bølgefunktion i Psi(x,0): ---
% % Gaussisk bølgefunktion:
a = 1/10;
x0 = -10;
Psi0 = exp(-a*(x-x0).^2);

% Moduler bølgefunktionen med en eksponential-faktor, så den får en
% middelimpuls forskellig fra nul:
k0 = 2;
Psi0 = Psi0.*exp(1i*k0*x);

% --- Normer bølgefunktionen: ---
Norm = sum(abs(Psi0).^2)*dx;  % Dette er bølgefunktionens norm-kvadrat.
Psi0 = Psi0/sqrt(Norm);  % Nu er Psi normeret.


% --- Beregn spredningen af x for denne bølgefunktion: ---
Meanx = sum(x.*abs(Psi0).^2)*dx;
MeanxSqr = sum(x.^2.*abs(Psi0).^2)*dx;
sigmax = sqrt(MeanxSqr - Meanx^2);

% Noter fornuftige plotte-grænser for det definerede x og for |Psi(x,t)|^2.
XLIM = [min(x) max(x)];
AbsPsi0 = abs(Psi0);
AbsSqrPsi0 = AbsPsi0.^2;
YLIM = [-1.4 1.4]*max(AbsPsi0);
YLIMSQR = [0 2]*max(AbsSqrPsi0);

% --- Definer en k-akse, hvorpå phi(k) skal tage sine værdier: -----
% Bemærk, hvis den angivne bølgefunktion er en
% minimum-usikkerheds-tilstand, så er sigmax*sigmap = 1/2 (da hbar = 1). Vi
% har også at p = k (igen da hbar = 1), så det er nok fornuftigt at
% definere phi(k) på følgende k-akse (ellers, lav din egen):
kmin = -5/sigmax + k0;
kmax = 5/sigmax + k0;
dk = 0.01/sigmax;
k = [kmin:dk:kmax];

% For hvert k skal vi beregne beta, B og F. Vi antaget A = 1 i ligning [2.137]:
beta = alpha./k;  % Ligning [2.135].
B = 1i*beta./(1-1i*beta); % [2.137].
F = 1./(1-1i*beta);       % [2.137].

% Beregn phi(k) fra ligning [2.103] - eller en generellisering af [2.103]:
phi = zeros(size(k)); % Alloker plads til phi(k).
for j = [1:length(k)]
    k_now = k(j);  % Den aktuelle værdi af k.
    B_now = B(j);
    F_now = F(j);
    %     psi_k_now = 1/sqrt(2*pi)*((exp(1i*k_now*x) + B_now*exp(-1i*k_now*x)).*heaviside(-x) ...   % Fører til en generelisering
    %                               + F_now*exp(1i*k_now*x).*heaviside(x));                         % af [2.103].
    psi_k_now = 1/sqrt(2*pi)*exp(1i*k_now*x);   % Svarer til "indmaden" af integralet i [2.103].
    % Prøv begge af de ovenstående metoder. Er der forskel? Hvorfor/hvorfor ikke?
    phi(j) = sum(Psi0.*conj(psi_k_now))*dx;  % Nu udføres integralet fra [2.103] (eller genereliseringen deraf).
end


% Angiv nu et antal tidspunkter, hvor bølgefunktionen ønskes plottet
% juster selv efter behov:
tstep = 0.2;
tmax = 2*abs(x0/k0); % Valgt således, at midten af den reflekterede puls slutter samme sted, som den indkommende puls startede.
t = [0:tstep:tmax];

% --- Gør figuren klar til plotning:
figure(1);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).

for j = [1:length(t)]
    for l = [1:length(x)]
        t_now = t(j); % Den aktuelle værdi af t.
        x_now = x(l); % Den aktuelle værdi af x.
        if x_now <= 0
            Psi(l) = sum(phi.*(exp(1i*k*x_now) + B.*exp(-1i*k*x_now)).*exp(-1i*k.^2/2*t_now))*dk; % [2.137]
        else
            Psi(l) = sum(phi.*F.*exp(1i*k*x_now).*exp(-1i*k.^2/2*t_now))*dk; % [2.137]
        end
        % Ovenstående svarer til [2.100], men genereliseret til
        % delta-funktionspotentialet.
    end
    Psi = Psi/sqrt(2*pi); % Vi manglede lige denne kvadratrod ovenfor.
    
    AbsSqrPsi = abs(Psi).^2; % Denne vil vi plotte.
    
    plot(x,AbsSqrPsi,'k-',x,AbsSqrPsi0,'k:',[0 0],YLIM,'k--');
    title(['t = ' sprintf('%.3f',t_now) '.']);
    axis([XLIM YLIMSQR]);
    
    xlabel('x');
    ylabel('|\Psi(x,t)|^2');
    drawnow;
end

% -----------------------------------------------------------------
% Beregn en reflektions og transmissionskoefficient for HELE PULSE:
% -----------------------------------------------------------------
idxm = find(x < 0);
idxp = find(x > 0);
Area_in = sum(AbsSqrPsi0(idxm))*dx;   % Integral under dashed curve.
Area_ref = sum(AbsSqrPsi(idxm))*dx;   % Integral under reflected pulse.
Area_trans = sum(AbsSqrPsi(idxp))*dx; % Integral under transmitted pulse.
% Bemærk, for at ovenstående virker, så skal beregningerne i figur 1 slutte
% med vel-separerede pulser, som ikke må have forladt det aktuelle
% plotteområde.


% Vi gentager nu proceduren fra figur 1, men plotter istedet real- og
% imaginær-del.

% --- Gør figuren klar til plotning:
figure(2);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).

for j = [1:length(t)]
    for l = [1:length(x)]
        t_now = t(j); % Den aktuelle værdi af t.
        x_now = x(l); % Den aktuelle værdi af x.
        if x_now <= 0
            Psi(l) = sum(phi.*(exp(1i*k*x_now) + B.*exp(-1i*k*x_now)).*exp(-1i*k.^2/2*t_now))*dk; % [2.137]
        else
            Psi(l) = sum(phi.*F.*exp(1i*k*x_now).*exp(-1i*k.^2/2*t_now))*dk; % [2.137].
        end
    end
    Psi = Psi/sqrt(2*pi);
    
    RePsi = real(Psi); % Denne vil vi plotte.
    ImPsi = imag(Psi); % Denne vil vi plotte.
    AbsPsi = abs(Psi); % Denne vil vi plotte.
    if j==1
        h = plot(x,RePsi,'b-',x,ImPsi,'r-',x,AbsPsi,'k-',x,AbsPsi0,'k:',[0 0],YLIM,'k--');
        title(['t = ' sprintf('%.3f',t_now) '.']);
        axis([XLIM YLIM]);
        xlabel('x');
        ylabel('\Psi(x,t)');
        legend('Re','Im','Abs','location','SouthWest');
        drawnow;
    else
        set(h(1),'ydata',RePsi)
        set(h(2),'ydata',ImPsi)
        set(h(3),'ydata',AbsPsi)
        title(['t = ' sprintf('%.3f',t_now) '.']);
        drawnow
    end
end

% Plot R og T som funktion af k, og plot også værdierne for PULSERNE:
figure(3);
hold off;
plot(k,abs(phi).^2,'k:');
hold on;

R = beta.^2./(1+beta.^2);  % [2.138]
T = 1./(1+beta.^2);        % [2.139]
plot(k,R,'b-');
plot(k,T,'r-');

R_pulse = Area_ref/Area_in*ones(size(k));
T_pulse = Area_trans/Area_in*ones(size(k));
plot(k,R_pulse,'b--');
plot(k,T_pulse,'r--');

xlabel('k');
legend('|\phi(k)|^2','R(k)','T(k)','R(pulse)','T(pulse)','location','EastOutside');