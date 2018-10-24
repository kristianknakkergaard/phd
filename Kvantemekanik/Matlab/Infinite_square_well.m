% ------------------------------------------------------------------
%  Denne Matlab-kode skal illustrere forskellige l�sninger til den
%  uendelige br�nd, som beskrevet i kapitel 2.2 i Griffith's "Introduction
%  to Quantum Mechanics".
%
%  Bem�rk, vi har sat hbar = m = 1 i de aktuelle formler.
% ------------------------------------------------------------------
close all
clear all
FontSize = 14;
set(0,'defaulttextinterpreter','latex');

% --- Definer f�rst en x-akse: ---
a = 1;
x = [0:0.001:1]*a;
dx = 0.001*a;

% --- Definer s� en start-b�lgefunktion Psi(x,0): ---
% Husk, den skal v�re nul i x=0 og x=a:
Psi0 = zeros(size(x)); % Nul over alt.
idx = find((x >= 0) & (x <= a)); % Index til omr�det inde i br�nden.

% Kommenter og udkommenter nedenfor efter behov, eller indtast din egen
% startb�lgefunktion Psi0:

% Eksempel 2.2
% Psi0(idx) = x(idx).*(a-x(idx)); % Bem�rk: ikke normeret.

% % Modificeret Eksempel 2.2, forceret til nul i x=a/2.
% Psi0(idx) = x(idx).*(a-x(idx)).*(a/2-x(idx)); % Bem�rk: ikke normeret.

% % Eksempel 2.2, men kvadreret.
% Psi0(idx) = x(idx).^2.*(a-x(idx)).^2; % Bem�rk: ikke normeret.

% % Trekant:
 idx1 = find((x >= 0) & (x <= a/2));
 idx2 = find((x > a/2) & (x <= a));
 Psi0(idx1) = x(idx1);
 Psi0(idx2) = a-x(idx2);  % Bem�rk, ikke normaliseret.



% Pr�v evt. at modulere b�lgefunktionen med en eksponential-faktor
% (s�t k0 til noget andet end nul i linjen nedenfor):
k0 = 0;
Psi0 = Psi0.*exp(1i*k0*x);

% --- Normer b�lgefunktionen: ---
Norm = sum(abs(Psi0).^2)*dx;  % Dette er b�lgefunktionens norm-kvadrat.
Psi0 = Psi0/sqrt(Norm);  % Nu er Psi normeret.


% Noter fornuftige plotte-gr�nser for det definerede x og for |Psi(x,t)|^2.
XLIM = [0 a];
AbsPsi0 = abs(Psi0);
AbsSqrPsi0 = AbsPsi0.^2;
YLIM = [-1.4 1.4]*max(AbsPsi0);
YLIMSQR = [0 1.4]*max(AbsSqrPsi0);


% Angiv nu et antal tidspunkter, hvor b�lgefunktionen �nskes plottet
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

Sum_cn_square = 0; % Hold styr p� koefficienterne - skal give et til sidst.
while Sum_cn_square < (1 - 1e-7)
    psin = zeros(size(x)); % Nul over alt.
    psin(idx) = sqrt(2/a)*sin(n*pi*x(idx)/a); % Nu er psi_n(x) defineret.
    cn = sum(psin.*Psi0)*dx; % Dette er ligning [2.37] eller [2.34].
    Sum_cn_square = Sum_cn_square + abs(cn)^2; % Opdater aktuel sum - skal give 1 jvf. [2.38].
    
    % Opdater vektorer:
    cn_vec = [cn_vec; cn]; % Inds�t nyeste cn.
    n_vec = [n_vec; n]; % Inds�t nyeste n.
    n = n + 1;
end



% --- G�r figuren klar til plotning:
figure(1);
set(gca,'FontSize',FontSize);


for j = [1:length(t)]
    Psi = zeros(size(x)); % Alloker plads til Psi(x,t) og s�t til nul.
    t_now = t(j); % Den aktuelle v�rdi af t i ligning [2.36].
    for k = [1:length(n_vec)]  % L�b igennem de aktuelle n-v�rdier.
        n = n_vec(k); % Aktuelle n-v�rdi.
        cn = cn_vec(k); % Aktuelle cn-v�rdi.
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


% --- G�r figuren klar til plotning:
figure(2);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).

for j = [1:length(t)]
    Psi = zeros(size(x)); % Alloker plads til Psi(x,t) og s�t til nul.
    t_now = t(j); % Den aktuelle v�rdi af t i ligning [2.36].
    for k = [1:length(n_vec)]  % L�b igennem de aktuelle n-v�rdier.
        n = n_vec(k); % Aktuelle n-v�rdi.
        cn = cn_vec(k); % Aktuelle cn-v�rdi.
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

