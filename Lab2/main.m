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


%% Problem 5: Test av normalitet
%{
    - Visuell undersökning ger:
        - Moderns ålder är någorlunda N-fördelade men något höger skev.
        - Moderns längd ser ut att vara N-fördelad, med viss koncentration
          kring 160 cm och 170 cm.
        - Moderns vikt avviker från N-fördelningen vid övre extremen, höger skev.
        - Barnets vikt avviker från N-fördelningen vid lägre extremen, vänster skev.
    
    - Jarque-Beras test vid nivå 5% ger:
        - Att moderns ålder är N-fördelad förkastas 
        - Att moderns längd är N-fördelad förkastas inte
        - Att moderns vikt är N-fördelad förkastas 
        - Att barnets vikt är N-fördelad förkastas
%}
clc; clear variables; clf; close all;
load birth.dat

age = birth(:, 4);          % Moderns ålder i år
height = birth(:, 16);      % Moderns längd i cm
weight = birth(:, 15);      % Moderns vikt i kg
childWeight = birth(:, 3);  % Barnets vikt i g

% Plotta jämförelse med N-fördelning
figure(1)
subplot(2,2,1)
normplot(age), title('Moderns ålder [år]')
subplot(2,2,2)
normplot(height), title("Moderns längd [cm]")
subplot(2,2,3)
normplot(weight), title("Moderns vikt [kg]")
subplot(2,2,4)
normplot(childWeight), title("Barnets vikt [g]")

% Använd Jarque-Beras test för att avgöra om N-fördelat
% H_0: att datan är N-fördelad
% H_1: datan är inte N-fördelad
% 1 -> förkasta H_0 
% 0 -> förkasta inte H_0
alpha = 0.05;
ageJB = jbtest(age, alpha);
heightJB = jbtest(height, alpha);
weightJB = jbtest(weight, alpha);
childWeightJB = jbtest(childWeight, alpha);

% Plotta täthetsfunktionerna
% Bara för att se hur de ser ut
figure(2)
subplot(2,2,1)
ksdensity(age), title('Moderns ålder [år]')
subplot(2,2,2)
ksdensity(height), title("Moderns längd [cm]")
subplot(2,2,3)
ksdensity(weight), title("Moderns vikt [kg]")
subplot(2,2,4)
ksdensity(childWeight), title("Barnets vikt [g]")

% Skriv ut resultat av Jarque-Beras testet
fprintf('Följande gäller för signifikansnivån %d%%:', alpha*100)
fprintf('\nAtt moderns ålder är N-fördelad förkastas '), if ageJB==0, fprintf("inte"), end
fprintf('\nAtt moderns längd är N-fördelad förkastas '), if heightJB==0, fprintf("inte"), end
fprintf('\nAtt moderns vikt är N-fördelad förkastas '), if weightJB==0, fprintf("inte"), end
fprintf('\nAtt barnets vikt är N-fördelad förkastas '), if childWeightJB==0, fprintf("inte"), end

%% Problem 6: Enkel linjär regression - Moores lag
%{
    - Undersökning av Moores lag för transistorers ytdensitet.
    - Residualen ser ut att vara N-fördelad, dock något vänster skev.
    
    Info om linjär regression:
    - https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
%}
clc; clear variables; clf; close all; format long g;
load moore.dat

x = moore(:, 1); % år
y = moore(:, 2); % transistorer/ytenhet

w = log(y); % logaritmerar för att gå från exponentiellt till linjärt

X = [ones(length(x),1), x]; % Bilda matris för att använda i regress

% Detta är samma som X\w tror jag
% X\w är MK-skattning och jag tror regress gör samma i det här fallet
% STATS är en vektor där första elementet är R^2
[B_hat,~,~,~,STATS] = regress(w, X); % skattning av [beta_0; beta_1]

w_hat = X*B_hat; % skattad modell, w_hat=log(y_hat)=X*B_hat

% Plotta skattningen och riktiga datan
figure(1)
plot(x, w)
hold on
plot(x, w_hat)
legend('Riktiga datan', 'Skattning', 'Location', 'best')
xlabel('År'), ylabel('log(transistorer/ytenhet)')

% Plotta residualen (observerat värde - modellens värde)
% Ser ut som en något vänster skev N-fördelning
figure(2)
res = w-w_hat;
subplot(2,1,1), normplot(res)
subplot(2,1,2), hist(res)

% Bestämning av R^2 (Coefficient of determination)
% Varierar mellan 0-1 och ju högre värdet är desto bättre är modellen
% (Beräknas enligt: Rsq = 1 - sum((w - w_hat).^2)/sum((w - mean(w)).^2))
Rsq = STATS(1);

% Funktion för antal transistorer/ytenhet för något år enligt modellen
yYear = @(year) exp(B_hat(1) + B_hat(2)*year);
year = 2025;
fprintf('Prediktion för antalet transistorer år %d: %.0f st/ytenhet \n', year, yYear(2025))
fprintf('R^2 är %.5f \n', Rsq)





