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
