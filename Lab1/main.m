%% Problem 1: Simulering av exponentialfördelade slumptal
clear variables; clc; clf;
%{ 
    Genererar N stycken Exp(1/10)-fördelade slumptal
    Ritar upp ett histogram av slumptalen
    Plottar den sanna täthetsfunktionen som jämförelse med histogrammet

    Det är viss skillnad mellan histogrammet och täthetsfunktionen
    eftersom histogrammet inte är helt exakt beräknat
%}
mu = 10;
N = 1e4;
y = exprnd(mu, N, 1);           % Genererar N exp-slumptal
hist_density(y);                % Skapar ett normaliserat histogram
t = linspace(0, 100, N/10);     % Vektor med N/10 punkter
hold on
plot(t, exppdf(t, mu), 'r')     % exppdf är täthetsfunktionen
hold off

%% Problem 2: Stora talens lag
clear variables; clc; clf;
%{ 
    Plottar linje av rikiga väntevärdet
    Beräknar medelvärdet av n antal exponetialfördelande stokastiska
    variabler där n börjar vid 1 och går till 500 och plottar den punkten
    Man ser då att när n blir större så närmar sig medelvärdet riktiga
    väntevärdet
%}
mu = 0.5;
M = 500;
X = exprnd(mu, M, 1);
plot(ones(M, 1)*mu, 'r-.')
hold on
for k = 1:M
    plot(k, mean(X(1:k)), '.')
%     % Det här funkar men suger för hastighet
%     if k == 1         
%         % Lägg bara till legend först loopen och ta inte med data1...dataM
%         legend('Sant \mu', 'Skattning av \mu', 'AutoUpdate', 'off')
%     end
    xlabel(num2str(k)), pause(0.001)
end
legend('Sant \mu', 'Skattning av \mu')

%% Problem 3: Väntevarde av exp.fördelad stokastisk variabel
clear variables; clc;
%{ 
    Sätt väntevärdet mu till vad som helst och variera N
    Medelvärdet av N antal exponetialfördelande stokastiska
    variabler blir då relativt nära mu om N är stort och inte så nära om N
    är litet
%}
mu = 15;                
N = 1e4;
y = exprnd(mu, N, 1);
mean(y)

%% Problem 4: Monte Carlo-skattning av talet pi
clear variables; clc; clf;
%{
    Generera vektorerna U och V med slumptal i intervallet (-1, 1) och
    plotta dessa punkter
    Plotta enhetscirkeln
    Generera vektorn Z där varje element får värdet 1 om det ligger
    innanför enhetscirkeln och 0 om det ligger utanför, beräkna sedan
    medelvärdet och gångra med 4 för att får en skattning av pi
    
    Öka N för en mer noggrann skattning
%}

N = 1e3;
U = 2*rand(1,N)-1;          % Genererar U(-1,1)-fördelade slumptal
V = 2*rand(1,N)-1;
plot(U,V,'o'), hold on      % Plottar de genererade punkterna
X = -1:0.01:1;
plot(X, sqrt(1-X.^2), 'r')  % Plottar enhetscirkeln
plot(X, -sqrt(1-X.^2), 'r')
Z = (sqrt(U.^2+V.^2) <= 1); % Beräknar närmevärde på pi
pi = 4*mean(Z)

%% Problem 5: Beräkning av sannolikheter
clear variables; clc;
%{
    Allmänt: 
        XXXcdf är fördelningsfunktionen [ F_X(x), P(X <= x) ]
        XXXpdf är täthetsfunktionen [ f_X(x) ]          
%}

P_X1_1 = binocdf(3, 10, 0.3)
P_X1_2 = 1 - binocdf(7, 10, 0.3)
P_x1_3 = binocdf(4, 10, 0.3) - binocdf(3, 10, 0.3)

P_X2_1 = normcdf(3, 5, 3)
P_X2_2 = 1 - normcdf(7, 5, 3)
P_x2_3 = normcdf(4, 5, 3) - normcdf(3, 5, 3)

lambda_exp = 7;
mu_exp = 1/lambda_exp; % Väntevärdet för exp. fördelning
P_X3_1 = expcdf(3, mu_exp)
P_X3_2 = 1 - expcdf(7, mu_exp)
P_x3_3 = expcdf(4, mu_exp) - expcdf(3, mu_exp)

%% Problem 6: Visualisering av sannolikhetsfördelningar
% Problem 6: Täthetsfunktion för normalfördelning
clear variables; clc; clf;
mu = 0;
sigma = 1;
dx = 0.01;
x = -10:dx:10; 
y = normpdf(x, mu, sigma);
figure(1)
plot(x,y)

% Problem 6: Täthetsfunktion för gammafördelning
dx = 0.01;
x = 0:dx:10; 
y = gampdf(x,1,2);
figure(2)
plot(x,y), hold on
z = gampdf(x,5,1);
plot(x,z,'r')

% Problem 6: Fördelningsfunktion för gammafördelning
dx = 0.01;
x = 0:dx:10;
y = gamcdf(x,1,2);
figure(3)
plot(x,y), hold on
z = gamcdf(x,5,1);
plot(x,z,'r')






