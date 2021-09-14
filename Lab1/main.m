%% Problem 1: Simulering av exponentialfordelade slumptal
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