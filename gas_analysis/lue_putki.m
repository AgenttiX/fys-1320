% FYS-1320 Gas Analysis
% Copyright Alpi Tolvanen and Mika Mäki, 2016
% This code has been developed on Octave and
% it has not been tested on Matlab.

% Mean = keskiarvot sisältävä matriisi
% kansio: polku, joka loppuu /-merkkiin
% Tiedostonimet: cell array, joka sisältää string-muotoiset tiedostonimet

% KKV = keskiarvon keskivirhe
function [Mean, TempR_mean, KKV_max, KKV_max_rel, TempR_KKV_max, TempR_KKV_max_rel] = lue_putki(kansio, Tiedostonimet)

% Käytetään vain ensimmäisen minuutin mittausdataa
pisteita = 60;
suotimia = length(Tiedostonimet(:,1));

% Muodostetaan tensorit (moniulotteiset matriisit) datalle
Data = zeros(suotimia , 2, pisteita);
% Data2 = zeros(suotimia , 2, pisteita);
T_Data = zeros(suotimia , 2, pisteita);

for y = 1:suotimia;
    % Luetaan data matriiseiksi
    mittaus = pilkut(load(strcat(kansio, Tiedostonimet{y,1} )));
    referenssi = pilkut(load(strcat(kansio, Tiedostonimet{y,2} )));
    
    % Otetaan matriiseista peak-to-peak-vektori ja tallennetaan se tensoriin
    Data(y,1,:) = mittaus(:, 2);
    Data(y,2,:) = referenssi(:, 2);
    
    # Data2(y,1,:) = mittaus(:,3);
    # Data2(y,2,:) = referenssi(:,3);
    
    T_Data(y,1,:) = mittaus(:,4);
    T_Data(y,2,:) = referenssi(:,4);
end

Mean = zeros(suotimia, 2);
% Mean2 = zeros(suotimia, 2);
TempR_mean = zeros(suotimia, 2);
KKV = zeros(suotimia, 2);

for y = 1:suotimia;
    for x = 1:2;
        % 1. kanava
        Mean(y,x) = mean(Data(y,x,:));
        
        % 2. kanava
        % Mean2(y,x) = mean(Data2(y,x,:));
        
        % Termistorin resistanssi
        TempR_mean(y,x) = mean(T_Data(y,x,:));
        
        % Keskiarvon keskivirhe
        KKV(y,x) = std(Data(y,x,:)) ./ sqrt(pisteita);
        
        % Termistorin keskiarvon keskivirhe
        TempR_KKV(y,x) = std(T_Data(y,x,:)) ./ sqrt(pisteita);
    end
end

% Lasketaan yksittäisiä lukuarvoja
KKV_max = max(max(KKV));
TempR_KKV_max = max(max(TempR_KKV));

KKV_rel = KKV ./ Mean;
KKV_max_rel = max(max(KKV_rel));

TempR_KKV_rel = TempR_KKV ./ TempR_mean;
TempR_KKV_max_rel = max(max(TempR_KKV_rel));

end
