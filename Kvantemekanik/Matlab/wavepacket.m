% ------------------------------------------------------------------
%  Denne Matlab-kode skal illustrere forskellige løsninger til den fri
%  partikel, som beskrevet i kapitel 2.4 i Griffith's "Introduction to
%  Quantum Mechanics".
%
%  Bemærk, vi har sat hbar = m = 1 i de aktuelle formler. Så vil
%  bølgetallet k og impulsen p have samme enheder.
% ------------------------------------------------------------------
close all
clear all
FontSize = 14;

% Vis analytisk resultat for psi(x,t)? 1 for ja, 0 for nej.
% Stemmer kun mening, hvis en gaussisk start Psi0 er valgt og k0 = 0.
ShowGaussianAnalytical = 0;

% --- Definer først en x-akse: ---
xmin = -10;
xmax = 10;
dx = 0.05;
x = [xmin:dx:xmax];

% --- Definer så en start-bølgefunktion Psi(x,0): ---
% Kasse-bølgefunktion:
a = 1;
Psi0 = (abs(x) <= a);  % Nu er Psi = 1 i intervallet -a <= x <= a.

% % % Gaussisk bølgefunktion:
% a = 1;
% Psi0 = exp(-a*x.^2);

% til opg. 2.21:
% a = 1; 
% Psi0 = sqrt(a).* exp(-a.*abs(x));

% En halv periode af en cosinus-funktion:
% a = 3;
% idx = find((x >= -a/2) & (x <= a/2));
% Psi0 = zeros(size(x));
% Psi0(idx) = sqrt(2/a)*cos(pi*x(idx)/a);

% Prøv evt. at modulere bølgefunktionen med en eksponential-faktor
% (sæt k0 til noget andet end nul i linjen nedenfor):
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
kmin = -10/sigmax;
kmax = 10/sigmax;
dk = 0.05/sigmax;
k = [kmin:dk:kmax];

% --- Beregn phi(k) fra ligning [2.103]:
phi = zeros(size(k)); % Alloker plads til phi(k).
for j = [1:length(k)]
    k_now = k(j);  % Den aktuelle værdi af k.
    phi(j) = 1/sqrt(2*pi)*sum(Psi0.*exp(-1i*k_now.*x))*dx;  % Dette er [2.103].
end

% Beregn momenter af p:
Meank = sum(abs(phi).^2.*k)*dk;  % Middelværdi.
MeankSqr = sum(abs(phi).^2.*k.^2)*dk;  
sigmak = sqrt(MeankSqr - Meank^2);  % Spredning.




% Angiv nu et antal tidspunkter, hvor bølgefunktionen ønskes plottet
% juster selv efter behov:
tstep = 0.02;
tmax = 1;
%t = [-tmax:tstep:tmax];
t = [0:tstep:tmax 0];

% --- Gør figuren klar til plotning af sandsynlighedsfordelingen: -------
figure(1);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).

for j = [1:length(t)]
    for l = [1:length(x)]
        t_now = t(j); % Den aktuelle værdi af t i ligning [2.100].
        x_now = x(l); % Den aktuelle værdi af x i ligning [2.100].
        Psi(l) = sum(phi.*exp(1i*(k*x_now - k.^2/2*t_now))*dk); % Dette er næsten [2.100].
    end
    Psi = Psi/sqrt(2*pi); % Så, nu har vi beregnet [2.100].

    AbsSqrPsi = abs(Psi).^2; % Denne vil vi plotte.

    plot(x,AbsSqrPsi,'k-',x,AbsSqrPsi0,'k:');
    
    Meanx_now = sum(AbsSqrPsi.*x)*dx;
    MeanxSqr_now = sum(AbsSqrPsi.*x.^2)*dx;
    sigmax_now = sqrt(MeanxSqr_now - Meanx_now.^2);
    Product = sigmax_now*sigmak;
    
    % Husk, k = p når hbar = 1.
    title(['t = ' sprintf('%.3f',t_now) ', \sigma_x = ' sprintf('%.3f',sigmax_now) ...
           ', \sigma_p = ' sprintf('%.3f',sigmak) ', \sigma_x*\sigma_k = ' sprintf('%.3f',Product) '.' ]);
       
    axis([XLIM YLIMSQR]);

    xlabel('x');
    ylabel('|\Psi(x,t)|^2');
    drawnow;
end


% --- Gør figuren klar til plotning af selve bølgefunktionen ----------
figure(2);
set(gca,'FontSize',FontSize);
Psi = zeros(size(x)); % Alloker plads til Psi(x,t).
Psi_theory = zeros(size(x)); % Alloker plads til den analytiske Psi(x,t) [problem 2.22].

for j = [1:length(t)]
    for l = [1:length(x)]
        t_now = t(j); % Den aktuelle værdi af t i ligning [2.100].
        x_now = x(l); % Den aktuelle værdi af x i ligning [2.100].
        Psi(l) = sum(phi.*exp(1i*(k*x_now - k.^2/2*t_now))*dk); % Dette er næsten [2.100].
    end
    Psi = Psi/sqrt(2*pi); % Så, nu har vi beregnet [2.100].
    
    RePsi = real(Psi); % Denne vil vi plotte.
    ImPsi = imag(Psi); % Denne vil vi plotte.
    AbsPsi = abs(Psi); % Denne vil vi plotte.
    [maxAbsPsi,posmax] = max(AbsPsi); % Her er sandsynlighedstætheden maximal.
    x0 = x(posmax);
    
    if ShowGaussianAnalytical
        Psi_theory = (2*a/pi)^(1/4)*exp(-a*x.^2/(1+2*1i*a*t_now))/sqrt(1+(2*1i*a*t_now));
        RePsi_theory = real(Psi_theory); % Denne vil vi plotte.
        ImPsi_theory = imag(Psi_theory); % Denne vil vi plotte.
        plot(x,RePsi_theory,'c.',x,ImPsi_theory,'g.', ...
           x,RePsi,'b-',x,ImPsi,'r-',x,AbsPsi,'k-',x,AbsPsi0,'k:',x0,0,'kx'); 
        legend('Re_t','Im_t','Re','Im','Abs','location','SouthWest');
    else
        plot(x,RePsi,'b-',x,ImPsi,'r-',x,AbsPsi,'k-',x,AbsPsi0,'k:',x0,0,'kx'); 
        legend('Re','Im','Abs','location','SouthWest');
    end
 
    title(['t = ' sprintf('%.3f',t_now) '.']);
    axis([XLIM YLIM]);

    xlabel('x');
    ylabel('\Psi(x,t)');
    drawnow;
end


% ----- Sammenlign phi(k) med det analytisk udtryk - husk at udkommentere
% det korrekte udtryk. Nedenstående udtryk antager en reel
% startbølgefunktion Psi(x,0).

phi_th = 1/sqrt(pi*a)*sin((k-k0)*a)./(k-k0);  % Gælder for en "kasse-bølgefunktion".
% phi_th = 1/(2*pi*a)^(1/4)*exp(-(k-k0).^2/4/a); % Gælder for en gaussisk bølgefunktion.
% phi_th = sqrt(4*a/pi^3)*cos((k-k0)*a/2)./(1 - ((k-k0)*a/pi).^2); % Gælder for "En halv periode af en cosinus-funktion"

figure(3);
hold off;
plot(k,real(phi),'bo');
hold on;
plot(k,imag(phi),'ro');
plot(k,real(phi_th),'c-', 'LineWidth',2);
plot(k,imag(phi_th),'k-', 'LineWidth',2);

legend('Re[\phi] nummerisk','Im[\phi] nummerisk','Re[\phi] analytisk','Im[\phi] analytisk');

