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


%% Problem 4: Fördelningar av givna data
%{
    - Plottar lådagram och uppskattade täthetsfunktioner för en datamängd
    - Undersöker hur rökning, motion och alkohol påverkar födelsevikt
    - För rökare kan man se en ganska tydlig minskning i födelsevikt
    - Lite motion verkar också i viss mån minska födelsevikten
    - Alkohol verkar inte ha någon tydlig påverkan på födelsevikt
%}
clc; clear variables; clf; close all;
load birth.dat

% ---- RÖKARE ----
% Kolla om mamman röker eller ej och plocka ut barnets vikt i x resp y
% (kolonn 20 är rökning och kolonn 3 är barnets vikt)
xSmoke = birth(birth(:, 20) < 3, 3);     % Barns vikt för mammor som INTE röker
ySmoke = birth(birth(:, 20) == 3, 3);    % Barns vikt för mammor som röker

figure('Name','Rökare')
% Plottar lådagram av datan
% Visar median, kvartiler och max min värden
% Notera att max min som stardard inte räknar med för stora avvikelser
% (värden som diffar mer än 1.5ggr avståndet mellan kvartilerna räknas inte med)
% (alltså totala variationsbredden som visas är 1.5+1+1.5 ggr av boxen)
% (För att ta med all data sätt 'whisker' till inf och då måste axis utökas)
subplot(2,2,1), boxplot(xSmoke), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Icke-rökare')
subplot(2,2,2), boxplot(ySmoke), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Rökare')

% Plottar uppskattningar av täthetsfunkterna för datan (ksdensity)
subplot(2,2,3:4) 
[fx, tx] = ksdensity(xSmoke);   % uppskattad täthetsfunktion för icke-rökare
plot(tx, fx, 'b')               % plottar denna separat för att vara konsekvent
hold on
[fy, ty] = ksdensity(ySmoke);   % uppskattad täthetsfunktion för rökare
plot(ty, fy, 'r')
legend('Icke-rökare', 'Rökare')
xlabel('Vikt [g]')
hold off

% ---- MOTION ----
% OBS! Det står fel i birth.txt, mycket motion är inte 2 utan 3 i datan
xExercise = birth(birth(:, 25) == 3, 3); % Mycket motion
yExercise = birth(birth(:, 25) == 1, 3); % Lite motion

figure('Name','Motion')
subplot(2,2,1), boxplot(xExercise), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Mycket motion')
subplot(2,2,2), boxplot(yExercise), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Lite motion')

subplot(2,2,3:4) 
[fx, tx] = ksdensity(xExercise); % uppskattad täthetsfunktion för mycket motion
plot(tx, fx, 'b')           
hold on
[fy, ty] = ksdensity(yExercise); % uppskattad täthetsfunktion för lite motion
plot(ty, fy, 'r')
legend('Mycket motion', 'Lite motion')
xlabel('Vikt [g]')
hold off

% ---- ALKOHOL ----
xAlcohol = birth(birth(:, 26) < 2, 3);  % Dricker inte alkohol
yAlcohol = birth(birth(:, 26) == 2, 3); % Dricker alkohol

figure('Name','Alkohol')
subplot(2,2,1), boxplot(xAlcohol), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Dricker inte alkohol')
subplot(2,2,2), boxplot(yAlcohol), 
axis([0 2 500 5000]), ylabel('Vikt [g]'), title('Dricker alkohol')

subplot(2,2,3:4) 
[fx, tx] = ksdensity(xAlcohol); % uppskattad täthetsfunktion för dricker inte
plot(tx, fx, 'b')           
hold on
[fy, ty] = ksdensity(yAlcohol); % uppskattad täthetsfunktion för dricker
plot(ty, fy, 'r')
legend('Ej alkohol', 'Alkohol')
xlabel('Vikt [g]')
hold off







