%% Problem 1: Simulering av konfidensintervall
%{
    - 95% av alla intervall förväntas innehålla mu 
      (i detta fall 95 stycken) vilket verkar stämma.
    - Det gröna vertikala strecket indikerar väntevärdet mu.
    - De horisontella strecken visar konfidensintervallen
      (de blå innehåller mu och de röda innehåller inte mu).
%}
clc; clear variables; clf;
% Parametrar:
n = 25;     % Antal mätningar
mu = 2;     % Väntevärdet
sigma = 1;  % Standardavvikelsen
alpha = 0.05;

% Simulerar n observationer för varje intervall
x = normrnd(mu, sigma, [n,100]); % n-by-100 matris med värden
% Skattar mu med medelvärdet
xbar = mean(x); % vektor med 100 medelvärden.
% Beräknar de undre och övre gränserna
undre = xbar - norminv(1-alpha/2)*sigma/sqrt(n);
ovre = xbar + norminv(1-alpha/2)*sigma/sqrt(n);

% Problem 1: Simulering av konfidensintervall (forts.)
% Ritar upp alla intervall
figure(1)
hold on
for k=1:100
    if ovre(k) < mu 
        % Rödmarkerar intervall som missar mu
        plot([undre(k) ovre(k)], [k k], 'r')
    elseif undre(k) > mu 
        % Rödmarkerar intervall som missar mu
        plot([undre(k) ovre(k)], [k k], 'r')
    else
        % Blåmarkerar intervall som innehåller mu
        plot([undre(k) ovre(k)], [k k], 'b')
    end
end
% b1 och b2 är bara till för att figuren ska se snygg ut
b1 = min(xbar - norminv(1 - alpha/2)*sigma/sqrt(n));
b2 = max(xbar + norminv(1 - alpha/2)*sigma/sqrt(n));
axis([b1 b2 0 101]) % Tar bort outnyttjat utrymme i figuren
% Ritar ut det sanna värdet
plot([mu mu], [0 101], 'g')
hold off


%% Problem 2: Maximum likelihood/Minsta kvadrat
%{
    - skattning av b med ML- och MK-skattning
    - båda skattningarna blir bra
    - täthetsfunktionen följer rayleighfördelningen bra
%}
clc; clear variables; clf;
M = 1e4;
b = 4;
x = raylrnd(b, [M, 1]);
hist_density(x, 40)
hold on
my_est_ml = sqrt(1/(2*length(x))*sum(x.^2));    % Räknad för hand
my_est_mk = sqrt(2/pi) * sum(x)/length(x);      % Räknad för hand
% Plotta skattningarna och b
plot(my_est_ml, 0, 'r*')
plot(my_est_mk, 0, 'g*')
plot(b, 0, 'ro')

% Plotta täthetsfunktionen
plot(0:0.1:10, raylpdf(0:0.1:10, my_est_ml), 'r')
hold off


%% Problem 3: Konfidensintervall för Rayleighfördelning
%{
    - MK-skatting och approximativt konfidensintervall
    - Konf. intervall gjort enligt paragraf 12.3
    - Medelfelet beräknat på papper (Förberedelseuppgift 2)
    - Täthetsfunktion passar fördelningen bra
%}
clc; clear variables; clf;
load wave_data.mat
subplot(2,1,1), plot(y(1:end))
subplot(2,1,2), hist_density(y)

alpha = 0.05;
n = length(y);

% MK-skattning från Problem 2
my_est = sqrt(2/pi) * sum(y)/n;
% Medelfel räknat på papper enligt Förberedelseuppgift 2
d = sqrt(2/pi * 1/n * (4-pi)/2 * my_est^2); 
% Beräknar de undre och övre gränserna (paragraf 12.3)
lower_bound = my_est - norminv(1-alpha/2)*d;
upper_bound = my_est + norminv(1-alpha/2)*d;

hold on
% Plotta skattningen och gränserna
plot(lower_bound, 0, 'g*')
plot(upper_bound, 0, 'g*')
plot(my_est, 0, 'ro')

% Plotta täthetsfunktionen
plot(0:0.1:6, raylpdf(0:0.1:6, my_est), 'r')
hold off

disp(['Approximativt ', num2str(100*(1-alpha)), '% konfidens intervall'])
disp(['my_est (MK): ', num2str(my_est)]) 
disp(['Undre gräns: ', num2str(lower_bound)]) 
disp(['Övre gräns: ', num2str(upper_bound)]) 







