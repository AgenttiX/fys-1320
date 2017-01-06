% Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

% Licensed with Creative Commons Attribution 4.0 International
% https://creativecommons.org/licenses/by/4.0/
% However, the license of the example code our work is based on is unclear and thereby so is the license
% for those parts of this code that are based on it

%9) ekstinkiton eli vaimenemisen laskeminen hiukkaskoon funktiona

N = 10000*10^5     % oletetaan että pitoisuus 10000/cm^3 -->10^9/m^3
L = 1         % oletetaan että "noin metrin pituinen putki" on metrin pitkä


% Käytetään pohjana ESIMERKKI_EKSTINKTIOTEHOKKUUS.m

% Esimerkkitiedosto ekstinktiotehokkuuksien laskentaan
% Joni Kalliokoski, 9.9.2014

m=1.33+0.001i; % veden taitekerroin
lambda=635; % valon aallonpituus (nm)
dp = logspace(1,4,1000); % hiukkaskokovektori, joissa ekstintiotehokuudet
                         % lasketaan (nm)

% Lasketaan Mie-parametrit:
Qext = zeros(size(dp)); %Alustetaan ekstintiovektori
for k = 1:numel(dp)
    result = Mie(m, pi*dp(k)/lambda); %Mie-parametrien laskenta
    Qext(k) = result(4); %Funktio palauttaa vektorin, jonka 4. alkio
    % on haluttu ekstinktiotehokkuus, muita arvoja ei tässä työssä tarvita.
end

sigma_ext=(pi*N.*(dp.*1e-9).^2.*Qext)/4;   %huomaa nanometrit

ext = 1- exp(-sigma_ext.*L);


%Piirretään ekstintiotehokkuus hiukkaskoon funktiona
figure(2)
loglog(dp,Qext,'-k')
xlabel('d_p (nm)')
ylabel('Q_{ext}')

%Piirretään ekstinktio hiukkaskoon funktiona
figure(1)
plot(dp,ext,'-k')
xlabel('d_p (nm)')
ylabel('ekstinktio')

