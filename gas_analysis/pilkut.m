% FYS-1320 Gas Analysis
% Copyright Alpi Tolvanen and Mika Mäki, 2016
% This code has been developed on Octave and
% it has not been tested on Matlab.

% Yhdistää tiedoston sarakkeet doubleiksi, sillä Octave
% lukee aluksi kokonaislukuosan ja desimaalit eri sarakkeisiin,
% koska tiedoston desimaalierotin on pilkku

% Tiedosto: 8-sarakkeinen matriisi
% Data: 4-sarakkeinen matriisi:
  % aika (min)
  % 1. kanavan peak-to-peak (mV)
  % 2. kanavan peak-to-peak (mV)
  % termistorin resistanssi (kOhm)

function Data=pilkut(Tiedosto);
    % y = rivinumero
    % y_v = rivi vektorina
    % for y=1:size(Tiedosto)[1];
    for y = 1:60;
        y_v=Tiedosto(y,:);
        
        % aika
        Data(y,1) = y_v(1) + y_v(2)*0.001;

        % 1. kanava
        Data(y,2) = y_v(3) + y_v(4)*0.001;
        
        % 2. kanava
        Data(y,3) = y_v(5) + y_v(6)*0.001;
        
        % termistori
        Data(y,4) = y_v(7) + y_v(8)*0.001;
    end;
end;
