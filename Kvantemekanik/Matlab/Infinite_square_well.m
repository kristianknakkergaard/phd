% ------------------------------------------------------------------
%  Denne Matlab-kode skal illustrere forskellige løsninger til den
%  uendelige brønd, som beskrevet i kapitel 2.2 i Griffith's "Introduction
%  to Quantum Mechanics".
%
%  Bemærk, vi har sat hbar = m = 1 i de aktuelle formler.
% ------------------------------------------------------------------
close all
clear all
FontSize = 14;
set(0,'defaulttextinterpreter','latex');

% --- Definer først en x-akse: ---
a = 1;
x = [0:0.001:1]*a;
dx = 0.001*a;

% --- Definer så en start-bølgefunktion Psi(x,0): ---
% Husk, den skal være nul i x=0 og x=a:
Psi0 = zeros(size(x)); % Nul over alt.
idx = find((x >= 0) & (x <= a)); % Index til området inde i brønden.

% Kommenter og udkommenter nedenfor efter behov, eller indtast din egen
% startbølgefunktion Psi0:

% Eksempel 2.2
% Psi0(idx) = x(idx).*(a-x(idx)); % Bemærk: ikke normeret.

% % Modificeret Eksempel 2.2, forceret til nul i x=a/2.
% Psi0(idx) = x(idx).*(a-x(idx)).*(a/2-x(idx)); % Bemærk: ikke normeret.

% % Eksempel 2.2, men kvadreret.
% Psi0(idx) = x(idx).^2.*(a-x(idx)).^2; % Bemærk: ikke normeret.

% % Trekant:
 idx1 = find((x >= 0) & (x <= a/2));
 idx2 = find((x > a/2) & (x <= a));
 Psi0(idx1) = x(idx1);
 Psi0(idx2) = a-x(idx2);  % Bemærk, ikke normaliseret.



% Prøv evt. at modulere bølgefunktionen med en eksponential-faktor
% (sæt k0 til noget andet end nul i linjen nedenfor):
k0 = 0;
Psi0 = Psi0.*exp(1i*k0*x);

% --- Normer bølgefunktionen: ---
Norm = sum(abs(Psi0).^2)*dx;  % Dette er bølgefunktionens norm-kvadrat.
Psi0 = Psi0/sqrt(Norm);  % Nu er Psi normeret.


% Noter fornuftige plotte-grænser for det definerede x og for |Psi(x,t)|^2.
XLIM = [0 a];
AbsPsi0 = abs(Psi0);
AbsSqrPsi0 = AbsPsi0.^2;
YLIM = [-1.4 1.4]*max(AbsPsi0);
YLIMSQR = [0 1.4]*max(AbsSqrPsi0);


% Angiv nu et antal tidspunkter, hvor bølgefunktionen ønskes plottet
% juster selv efter behov:

% Specifikt for eksempel 2.2:
E1 = pi^2/2/a^2;
E3 = pi^2*3^2/2/a^2;
T = 2*pi/(E3-E1);
t = [0:0.02:0.98]*T;

% Beregn cn'erne:
n = 1;
cn_vec = [];
n_vec = [];

Sum_cn_square = 0; % Hold styr på koefficienterne - skal give et til sidst.
while Sum_cn_square < (1 - 1e-7)
    psin = zeros(size(x)); % Nul over alt.
    psin(idx) = sqrt(2/a)*sin(n*pi*x(idx)/a); % Nu er psi_n(x) defineret.
    cn = sum(psin.*Psi0)*dx; % Dette er ligning [2.37] eller [2.34].
    Sum_cn_square = Sum_cn_square + abs(cn)^2; % Opdater aktuel sum - skal give 1 jvf. [2.38].
    
    % Opdater vektorer:
    cn_vec = [cn_vec; cn]; % Indsæt nyeste cn.
    n_vec = [n_vec; n]; % Indsæt nyeste n.
    n = n + 1;
end



% --- Gør figuren klar til plotning:
figure(1);
set(gca,'FontSize',FontSize);


for j = [1:length(t)]
    Psi = zeros(size(x)); % Alloker plads til Psi(x,t) og sæt til nul.
    t_now = t(j); % Den aktuelle værdi af t i ligning [2.36].
    for k = [1:length(n_vec)]  % Løb igennem de aktuelle n-værdier.
        n = n_vec(k); % Aktuelle n-værdi.
        cn = cn_vec(k); % Aktuelle cn-værdi.
        psin = zeros(size(x)); % Nul over alt.
        psin(idx) = sqrt(2/a)*sin(n*pi*x(idx)/a); % Nu er psi_n(x) defineret.
        Psi = Psi + cn*psin*exp(-1i*n^2*pi^2*t_now/2/a^2);  % Opdater med ligning [2.36].
    end

    AbsSqrPsi = abs(Psi).^2;
    
    plot(x,AbsSqrPsi,'k-',x,AbsSqrPsi0,'k:'); 
    title(['t = ' sprintf('%.3f',t_now) '.']);
    axis([XLIM YLIMSQR]);

    xlabel('$x$');
    ylabel('$|\Psi(x,t)|^2$');
    drawnow;
end


% Specifikt for eksempel 2.2:
E1 = pi^2/2/a^2;
T = 2*pi/E1;
t = [0:0.01:0.99]*T;


% --- Gør figuren klar til plotning:
figure(2);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).

for j = [1:length(t)]
    Psi = zeros(size(x)); % Alloker plads til Psi(x,t) og sæt til nul.
    t_now = t(j); % Den aktuelle værdi af t i ligning [2.36].
    for k = [1:length(n_vec)]  % Løb igennem de aktuelle n-værdier.
        n = n_vec(k); % Aktuelle n-værdi.
        cn = cn_vec(k); % Aktuelle cn-værdi.
        psin = zeros(size(x)); % Nul over alt.
        psin(idx) = sqrt(2/a)*sin(n*pi*x(idx)/a); % Nu er psi_n(x) defineret.
        Psi = Psi + cn*psin*exp(-1i*n^2*pi^2*t_now/2/a^2);  % Opdater med ligning [2.36].
    end

    RePsi = real(Psi); % Denne vil vi plotte.
    ImPsi = imag(Psi); % Denne vil vi plotte.
    AbsPsi = abs(Psi); % Denne vil vi plotte.

    plot(x,RePsi,'b-',x,ImPsi,'r-',x,AbsPsi,'k-',x,AbsPsi0,'k:'); 
    title(['t = ' sprintf('%.3f',t_now) '.']);
    axis([XLIM YLIM]);

    xlabel('$x$');
    ylabel('$\Psi(x,t)$');
    legend('Re','Im','Abs','location','SouthWest');
    drawnow;
end

