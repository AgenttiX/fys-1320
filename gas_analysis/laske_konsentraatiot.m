% FYS-1320 Gas Analysis
% Copyright Alpi Tolvanen and Mika Mäki, 2016
% This code has been developed on Octave and
% it has not been tested on Matlab.


% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.


% This file depends on .mat files provided in the FYS-1320 course materials.
% The rough outline of the matrix method was provided in the course
% example code, but our version has been rather heavily modified.

% The code is in Finnish and its quality decreases towards the end, but
% it works nonetheless.
% We are sorry for that but lack the motivation to fix it.


% Puhdistetaan työtila vanhoista muuttujista
clear;

% Määritellään vakiot
p = 1.013e5;        % paine
k = 1.3806505e-23;	% Boltzmannin vakio (J/K, MAOL)
T = 298;            % lämpotila (K) eli 25 C
L = 0.2;            % kaasunäytteen pituus (m)

% Tiedostokansiot
kansio1 = 'data/';
kansio2 = 'data2/';
vektorikansio = 'vektorit/';

% Ladataan pystyvektoreihin kaasujen absorbtiopinta-alat (sigma)
load(strcat(vektorikansio, 'kaasut_carbon.mat'), 'H2O', 'CO2', 'CO', 'CH4', 'C2H2', 'lambda');

% Huom! sigmat ovat yksikossa cm^2 -> muunnetaan yksikkoon m^2 kertoimella 1e-4
absorptio_matr = [H2O, CO2, CO, CH4, C2H2] * 1e-4; %vaakavektori

%Ladataan filtterien läpäisevyyskerroinvektorit
load(strcat(vektorikansio, '2710.mat'), 'NB2710');
load(strcat(vektorikansio, '3060.mat'), 'NB3060');
load(strcat(vektorikansio, '3220.mat'), 'NB3220');
load(strcat(vektorikansio, '3977.mat'), 'NB3977');
load(strcat(vektorikansio, '4270.mat'), 'NB4270');
load(strcat(vektorikansio, '4740.mat'), 'NB4740');
load(strcat(vektorikansio, '5756.mat'), 'NB5756');

% Suodinmatriisi
% (vaakavektori, joka sisältää pitkiä pystyvektoreita - kuten absorptio_matr)
suotimet=[NB2710, NB3060, NB3220, NB3977, NB4270, NB4740, NB5756]* 1e-4;

%Lisää tähän 0, 1 tai 2 suotimen järjestyslukua jotka laskenta ohittaa
poistettavat_suotimet=[];
% MIKÄLI LISÄTÄÄN, YKSITTÄISTEN KAASUJEN PIIRTO MENEE RIKKI 
% syystä että niitä kutsutaan järjestyluvulla 
% (matlabissa ei vektori-metodit oikein ole se vahvuus )

% Poistetaan halutut suotimet
for i = fliplr(sort(poistettavat_suotimet))
    suotimet(:,i)=[];        %Poistetaan esim. hiilidioksidi vektorin keskeltä
endfor

% Suodinten lukumäärä
s_lkm = size(suotimet)(2);



% Aallonpituusvälit lambda-vektorissa eivät ole vakioita, joten luodaan
% delta-lambda-vektori. Lisätään ensimmäiseksi väliksi toinen väli, jotta
% vektorin pituus säilyy.
% Välin pituus on aina positiivinen, joten otetaan itseisarvo.
dlambda_temp = diff(lambda);
dlambda = abs([dlambda_temp(1);dlambda_temp]);



% Luetaan mittaamamme arvot, matriisit ovat dimensiota 7*2,
% missä jokainen rivi sisältää irradianssin keskiarvon mitatulle suotimelle (:,1) 
% ja 3977-referenssille (:,2)

% Vertailuputki
M1_N2_nimet = {
'N2_2710-050.txt', 'N2_2710-050-reference.txt';
'N2_3060-060.txt', 'N2_3060-060_reference.txt';
'N2_3220-060.txt', 'N2_3220-060_reference.txt';
'N2_3977-065.txt', 'N2_3977-065.txt';
'N2_4270-070.txt', 'N2_4270-070_reference.txt';
'N2_4740-140.txt', 'N2_4740-140_reference.txt';
'N2_5756-080.txt', 'N2_5756-080_reference.txt';
};

[M1_N2_mean, M1_N2_TempR, M1_N2_KKV_max, M1_N2_KKV_max_rel, M1_N2_TempR_KKV_max M1_N2_TempR_KKV_max_rel] = lue_putki(kansio1, M1_N2_nimet);

% Näyteputki
M1_Sample_nimet = {
'Sample_2710-050.txt', 'Sample_2710-050_reference.txt';
'Sample_3060-060.txt', 'Sample_3060-060_reference.txt';
'Sample_3220-060.txt', 'Sample_3220-060_reference.txt';
'Sample_3977-060.txt', 'Sample_3977-060.txt';
'Sample_4270-070.txt', 'Sample_4270-070_reference.txt';
'Sample_4740-140_2.txt', 'Sample_4740-140_2_reference.txt';
'Sample_5756-080.txt', 'Sample_5756-080_reference.txt';
};

[M1_sample_mean, M1_sample_TempR, M1_sample_KKV_max, M1_sample_KKV_max_rel, M1_sample_TempR_KKV_max, M1_sample_TempR_KKV_max_rel] = lue_putki(kansio1, M1_Sample_nimet);


% ### 2. mittaus ###

% Vertailuputki
M2_N2_nimet = {
'M2_N2_2710-050.txt', 'M2_N2_3977-065.txt';
'M2_N2_3060-060.txt', 'M2_N2_3977-065.txt';
'M2_N2_3220-060.txt', 'M2_N2_3977-065.txt';
'M2_N2_3977-065.txt', 'M2_N2_3977-065.txt';
'M2_N2_4270-070.txt', 'M2_N2_3977-065.txt';
'M2_N2_4740-140.txt', 'M2_N2_3977-065.txt';
'M2_N2_5756-080.txt', 'M2_N2_3977-065.txt';
};

[M2_N2_mean, M2_N2_TempR, M2_N2_KKV_max, M2_N2_KKV_max_rel, M2_N2_TempR_KKV_max, M2_N2_TempR_KKV_max_rel] = lue_putki(kansio2, M2_N2_nimet);

% Hiukkasputki
M2_Hiu_nimet = {
'M2_Hiu_2710-050.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_3060-060.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_3220-060.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_3977-065.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_4270-070.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_4740-140.txt', 'M2_Hiu_3977-065.txt';
'M2_Hiu_5756-080.txt', 'M2_Hiu_3977-065.txt';
};

[M2_Hiu_mean, M2_Hiu_TempR, M2_Hiu_KKV_max, M2_Hiu_KKV_max_rel, M2_Hiu_TempR_KKV_max, M2_Hiu_TempR_KKV_max_rel] = lue_putki(kansio2, M2_Hiu_nimet);

% Ulkoputki
M2_out_nimet = {
'M2_out_2710-050.txt', 'M2_out_3977-065.txt';
'M2_out_3060-060.txt', 'M2_out_3977-065.txt';
'M2_out_3220-060.txt', 'M2_out_3977-065.txt';
'M2_out_3977-065.txt', 'M2_out_3977-065.txt';
'M2_out_4270-070.txt', 'M2_out_3977-065.txt';
'M2_out_4740-140.txt', 'M2_out_3977-065.txt';
'M2_out_5756-080.txt', 'M2_out_3977-065.txt';
};

[M2_out_mean, M2_out_TempR, M2_out_KKV_max, M2_out_KKV_max_rel, M2_out_TempR_KKV_max, M2_out_TempR_KKV_max_rel] = lue_putki(kansio2, M2_out_nimet);


disp(' ');
KKV_max_rel = [M1_N2_KKV_max_rel, M1_sample_KKV_max_rel, M2_N2_KKV_max_rel, M2_Hiu_KKV_max_rel, M2_out_KKV_max_rel]';
disp('Keskiarvon keskivirheen suhteelliset maksimit');
disp(KKV_max_rel);
disp('Maksimi:');
disp(max(KKV_max_rel));

disp(' ');
TempR_KKV_max_rel = [M1_N2_TempR_KKV_max_rel, M1_sample_TempR_KKV_max_rel, M2_N2_TempR_KKV_max_rel, M2_Hiu_TempR_KKV_max_rel, M2_out_TempR_KKV_max_rel]';
disp('Lämpötilan keskiarvon keskivirheen suhteelliset maksimit');
disp(TempR_KKV_max_rel);
disp('Maksimi:');
disp(max(TempR_KKV_max_rel));



disp(' ');
disp('1. mittauksen lisävirhearvio: referenssin keskiarvon keskivirhe mittausten välillä');

M1_N2_ref_KKV = std(M1_N2_mean(:,2)) ./ sqrt(s_lkm);
disp(M1_N2_ref_KKV);

M1_sample_ref_KKV = std(M1_sample_mean(:,2)) ./ sqrt(s_lkm);
disp(M1_sample_ref_KKV);

disp('Suhteellisena');
M1_N2_ref_KKV_rel = M1_N2_ref_KKV ./ mean(M1_sample_mean(:,2));
disp(M1_N2_ref_KKV_rel);

M1_sample_ref_KKV_rel = M1_sample_ref_KKV ./ mean(M1_sample_mean(:,2));
disp(M1_sample_ref_KKV_rel);

disp(' ');
disp('Typpiputken mittausten suhteellinen ero mittausten välillä (M2 / M1)');
disp( mean(mean(M2_N2_mean ./ M1_N2_mean)) );

% Seuraavat pätkät on kopioitu ja editoitu esimerkkikoodista


% Luodaan vektorit putkien irradiansseille kaikilla suotimilla
I_smp = M1_sample_mean(:,1);  % I_smp, ensimmäisen mittauksen tuntematon kaasunäyte
I0_M1 = M1_N2_mean(:,1);      % %I0_M1, tätä vastaavan typpi

I_Hiu = M2_Hiu_mean(:,1);     % I_Hiu, hiukkasen kiltahuoneen kaasunäyte
I_out = M2_out_mean(:,1);     % I_out, ulkoa otettu kaasunäyte
I0_M2 = M2_N2_mean(:,1);      % I0_M2,  toisen mittauskerran typpi

% Muistutuksena: sarake (:,1) tarkoittaa mitattuja arvoja, (:,2) mittauksen jälkeisiä referenssiarvoja



% Valonlähteen intensiteetin normalisointi
% Muutettu laskentatapaa 4.12 päivityksessä.
if(false) % Vanhatapa normalisoida
    I_smp = I_smp .* (M1_sample_mean(:,2)./M1_sample_mean(1,2));
    I0_M1 = I0_M1 .* (M1_N2_mean(:,2)./M1_N2_mean(1,2));

    I_Hiu = I_Hiu .* (M2_Hiu_mean(:,2)./M2_Hiu_mean(1,2));
    I_out = I_out .* (M2_out_mean(:,2)./M2_out_mean(1,2));
    I0_M2 = I0_M2 .* (M2_N2_mean(:,2)./M2_N2_mean(1,2));
else     % Uusitapa normalisoida
    I_smp = I_smp ./ (M1_sample_mean(:,2).*M1_N2_mean(4,2)./M1_sample_mean(4,2));
    I0_M1 = I0_M1 ./ (M1_N2_mean(:,2).*M1_N2_mean(4,2)./M1_N2_mean(4,2));

    I_Hiu = I_Hiu ./ (M2_Hiu_mean(:,2).*M2_N2_mean(4,2)./M2_Hiu_mean(4,2));
    I_out = I_out ./ (M2_out_mean(:,2).*M2_N2_mean(4,2)./M2_out_mean(4,2));
    I0_M2 = I0_M2 ./ (M2_N2_mean(:,2).*M2_N2_mean(4,2)./M2_N2_mean(4,2));
endif
% muutettu .* --> ./
% Lisätty termi "M1_N2_mean(4,2)", joka tekee siis sen että se palauttaa eri
% mittauskertojen putkien arvot samalle tasolle. Oikeastaan se ei tee mitää muuhun
% kuin mittauskertaan 1, kun käytimme vielä referenssejä hyväksemme. Joka tapauksessa
% Tämä pitäisi mielestäni olla virallinen mallikoodi koska nyt se on matemaattisesti tarkempi


% Käsitellään ylimääräisten suodinten poisto
for i = fliplr(sort(poistettavat_suotimet))
    I_smp(i)=[];         %Poistetaan haluttu suodin
    I0_M1(i)=[];

    I_Hiu(i)=[];
    I_out(i)=[];
    I0_M2(i)=[];
endfor



% Matriisiyhtälön ratkaiseminen

% ##### Lasketaan b #############
% s_lkm tarkoittaa suotimien lukumäärää, ja 'suotimet' ja 'dlamda' ovat korkeita matriiseja
% I ja I0 ovat pystyvektoreita, joiden pituus s_lkm, kuten transponoidun summalausekkeenkin
b_vec_smp=(1-I_smp./I0_M1) .* (sum( suotimet .*(dlambda*ones(1,s_lkm)) .*1e-2, 1))';  % Huom. transpoosi '

b_vec_Hiu=(1-I_Hiu./I0_M2) .* (sum( suotimet .*(dlambda*ones(1,s_lkm)) .*1e-2, 1))';

b_vec_out=(1-I_out./I0_M2) .* (sum( suotimet .*(dlambda*ones(1,s_lkm)) .*1e-2, 1))';




% ##### Lasketaan A #############

% Matriisin A sarakkeet
% Oikeellisuuden tarkistamiseksi voi verrata esimerkkikoodiin
sarake1_H2O = sum(absorptio_matr(:,1)*ones(1,s_lkm).*suotimet*1e-2.*(dlambda*ones(1,s_lkm)), 1)';  
sarake2_CO2 = sum(absorptio_matr(:,2)*ones(1,s_lkm).*suotimet*1e-2.*(dlambda*ones(1,s_lkm)), 1)';  
sarake3_CO  = sum(absorptio_matr(:,3)*ones(1,s_lkm).*suotimet*1e-2.*(dlambda*ones(1,s_lkm)), 1)';  
sarake4_CH4 = sum(absorptio_matr(:,4)*ones(1,s_lkm).*suotimet*1e-2.*(dlambda*ones(1,s_lkm)), 1)';  
sarake5_C2H2= sum(absorptio_matr(:,5)*ones(1,s_lkm).*suotimet*1e-2.*(dlambda*ones(1,s_lkm)), 1)';  

A = [sarake1_H2O, sarake2_CO2, sarake3_CO, sarake4_CH4, sarake5_C2H2].*p./(k*T)*L; %kerrottu jo aiemmin 1e-4.




% Konsentraatiot lasketaan pienimmän neliösumman menetelmällä
c_vec_smp = A\b_vec_smp;
c_vec_Hiu = A\b_vec_Hiu;
c_vec_out = A\b_vec_out;

% Muutetaan lukujen esitysmuotoa
format short g;

kaasut_ppm = [c_vec_smp, c_vec_Hiu ,c_vec_out]*1e6;

% Tulostetaan itse konsentraatiot ihan lopuksi.
#{
disp(' ');
disp('      Sample      Kiltahuone   Ulkoilma');
disp(kaasut_ppm);
disp('Järjestys: H20, CO2, CO, CH4, C2H2');
#}

% Tulostuvista riveistä toinen on siis hiilidioksidi, jolle lineaarinen approksimaatio ei toimi kovinkaan hyvin.
% Lisäksi tulosten oikeellisuuden epäileminen on perusteltua.





%-----------------------------------------------------------------
% Hiilidioksidin pitoisuuden laskeminen ja ensimmäisen kuvaajan valmistelu

% x-akseli
c_min = 0;
c_max = 1e-3;
arvauskonsentraatio=linspace(c_min,c_max,100);

% Lasketaan laskennallisesti  transmittanssi konsentraation funktiona (exp-käyrä)
transmittanssit_exp = sum(NB4270.*exp(-CO2.*(1e-4*p*L/(k*T))*arvauskonsentraatio).*dlambda*1e-2)/(sum(NB4270.*dlambda*1e-2));

% Lasketaan transmittanssi konsentraation funktiona (lineaarinen käyrä)
transmittanssit_lin = 1 - ((sum(NB4270.*CO2.*arvauskonsentraatio.*dlambda*1e-4)*(1e-2*p*L/(k*T))))/(sum(NB4270.*dlambda*1e-2));

% Mittaustuloksia kuvaavat vaakaviivat
mitattu_transmittanssi_smp = (I_smp(5)/(I0_M1(5)));
mitattu_transmittanssi_Hiu = (I_Hiu(5)/(I0_M2(5)));
mitattu_transmittanssi_out = (I_out(5)/(I0_M2(5)));

% Tässä lasketaan milloin ylemmät transmittanssit leikkaavat ekspotenttiaalisen käyrän
% Eli selvitetään konsentraatiot.
c_CO2_erot_nollasta_smp = transmittanssit_exp - mitattu_transmittanssi_smp;
c_CO2_erot_nollasta_Hiu = transmittanssit_exp - mitattu_transmittanssi_Hiu;
c_CO2_erot_nollasta_out = transmittanssit_exp - mitattu_transmittanssi_out;

% Etsitään nollakohdan indeksi. Tarkalleen ottaen graafinen määritys olisi tarkempi,
% koska se interpoloi alkioiden välistä, noh ero lienee jotain parin prosentin luokkaa
[i,indx_smp]=min(abs(c_CO2_erot_nollasta_smp));
[i,indx_Hiu]=min(abs(c_CO2_erot_nollasta_Hiu));
[i,indx_out]=min(abs(c_CO2_erot_nollasta_out));

% Kun tunnetaan indeksi niin tiedetään siitä sitten konsentraatio
c_CO2_numeerinen_smp = arvauskonsentraatio(indx_smp);
c_CO2_numeerinen_Hiu = arvauskonsentraatio(indx_Hiu);
c_CO2_numeerinen_out = arvauskonsentraatio(indx_out);

% Siirretty loppuun helppolukuisuuden vuoksi
#{
disp(' ');
disp('CO2:');
disp('      Sample      Kiltahuone  Ulkoilma');
disp([c_CO2_numeerinen_smp, c_CO2_numeerinen_Hiu, c_CO2_numeerinen_out]*1e6);
#}

% Vielä kuvaajan piirto

lw=2;        % linewidth
h1Fig=figure (1);
set(h1Fig, 'Position', [100 100 700 600]) %rikkoo asioita


plot(arvauskonsentraatio*1e6, transmittanssit_exp,  "linewidth", lw,
    arvauskonsentraatio(1,1:54)*1e6, transmittanssit_lin(1,1:54),  "linewidth", lw,
    [c_min,c_max]*1e6, [1,1].*mitattu_transmittanssi_smp,  "linewidth", lw,
    [c_min,c_max]*1e6, [1,1].*mitattu_transmittanssi_Hiu,  "linewidth", lw,
    [c_min,c_max]*1e6, [1,1].*mitattu_transmittanssi_out,  "linewidth", lw)
title("Transmittanssi CO2:n pitoisuuden funktiona NB4270-suotimen alueella",'FontSize', 20)
legendi=legend("Todellinen transmittanssi",
    "Transmittanssin lineaarinen approksimaatio",
    "Mitattu transmittanssi, laboratorion näyte",
    "Mitattu transmittanssi, Hiukkasen kiltahuone", 
    "Mitattu transmittanssi, ulkoilma");
set(legendi,'FontSize',18);
xlabel("Pitoisuus (ppm)",'FontSize', 20);
ylabel("Transmittanssi",'FontSize', 20);
axis([0,1000,0.6,1])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lisätty Sunnuntaina 4.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% TÄMÄ LISÄYS ON VAIN LISÄÄ VIRITTELYÄ, toisin sanoen laskentaa on muutettu hieman ja 
% lisää krääsää printaantuu näytölle. Ei ehkä ensisijaisesti HEAD:iksi, mutta
% ei siitä ei olisi haittaakaan.

% Tämä lisäys laskee matriisin virhetermejä (kohta-1) ja hiilimonoksidin yksittäisarvot (kohta-2)
% jälkimmäisestä piirretään myös kuvaaja

% Muutin myös itse normalisoinnin laskentatapaa, joka on nyt matemaattisesti tarkempi,
% joka vain ei toki korjaa mitään ongelmistamme.

% Ehkä vähän rumasti nyt siirsin itse konsentraatiomatriisin printtauksen loppuun
% mutta se johtuu siitä että koodiani en saa yhtään aiemmaksi, ja olisi siistimää
% jos matriisit printattaisiin ihan lopuksi.



% (kohta-1) PNS-virhetermien laskut vektorille b matriisista A:       (Ax=b)

% Sisältää pienimmän neliösumman ratkaisun suhteellisen residuaalin (itseisarvo), 
% eli suhteellisen virhetermin. Tätä voi ajatella että kun pns projisoi mittausdatasta
% tuloksiin, niin tämä kertoo kuinka paljon mittausdataa oli "muutettava" jotta 
% saatiin tulos jolla oli pienin virhe. Huomaa siis suhteellinen virhe, eli jaettu b:llä

% Niin elikkä siis:
% residuaali r = b - Ax  , missä x on pns ratkaisu, r residuaali (matikka 2 pruju sivu 9 3)
% virhe      norm(r)     , mikä on siis minimoitu
% ja omani:  abs(r./b)   , mikä oli täysin oma keksintöni katsoa minkä suotimien
%                          mittausarvoja "muutettiin" eniten jotta tulos paras

% Näihin menee tilaa sen verran että en tiedä kannattaako tehdä uutta diaa niitä varten,
% ainakin jos haluaa näyttää että me ollaan tehty jotain

format short

suhteellinen_virhe_smp = abs((b_vec_smp - A*c_vec_smp) ./ b_vec_smp);
suhteellinen_virhe_Hiu = abs((b_vec_Hiu - A*c_vec_Hiu) ./ b_vec_Hiu);
suhteellinen_virhe_out = abs((b_vec_out - A*c_vec_out) ./ b_vec_out);
virhe_kok_smp          = norm((b_vec_smp - A*c_vec_smp));
virhe_kok_Hiu          = norm((b_vec_Hiu - A*c_vec_Hiu));
virhe_kok_out          = norm((b_vec_out - A*c_vec_out));
disp(' ');
disp('Suotimille tehdyt suhteelliset korjauksen pns-menetelmässä');
disp('   Sample      Kiltahuone  Ulkoilma');
disp([suhteellinen_virhe_smp, suhteellinen_virhe_Hiu, suhteellinen_virhe_out]) % 7 suodinta
disp('absoluuttinen virhe b:ssä joka on minimoitu, (residuaalin normi)')
disp([virhe_kok_smp,virhe_kok_Hiu,virhe_kok_out]) % residuaalin normi

format short g;


% (kohta-2)
% Tämä laskee hiilimonoksidille saman kuin hiilidioksidille, eli tulos ilman linearisointia
% ja vielä kuvaajan piirto päälle. Kuvaajasta näkee tosi hyvin kuinka mitatut
% arvot ovat pyllyllään. Muuten koodi on kopioitu hiilidioksidin laskusta. 
% Kommentit poistettu turhaa tilaa vievinä.

if (length(poistettavat_suotimet) == 0) 
% iffin toinen pää 170 riviä myöhempänä, toisin sanoen yksittäisiä kuvaajia ei piirretä

CO_c_min = 0;
CO_c_max = 4e-3;
CO_arvauskonsentraatio=linspace(CO_c_min,CO_c_max,100);

CO_transmittanssit_exp = sum(NB4740.*exp(-CO.*(1e-4*p*L/(k*T))*CO_arvauskonsentraatio).*dlambda*1e-2)/(sum(NB4740.*dlambda*1e-2));
CO_transmittanssit_lin = 1 - ((sum(NB4740.*CO.*CO_arvauskonsentraatio.*dlambda*1e-4)*(1e-2*p*L/(k*T))))/(sum(NB4740.*dlambda*1e-2));
CO_mitattu_transmittanssi_smp = (I_smp(6)/(I0_M1(6)));
CO_mitattu_transmittanssi_Hiu = (I_Hiu(6)/(I0_M2(6)));
CO_mitattu_transmittanssi_out = (I_out(6)/(I0_M2(6)));
c_CO_erot_nollasta_smp = CO_transmittanssit_exp - CO_mitattu_transmittanssi_smp;
c_CO_erot_nollasta_Hiu = CO_transmittanssit_exp - CO_mitattu_transmittanssi_Hiu;
c_CO_erot_nollasta_out = CO_transmittanssit_exp - CO_mitattu_transmittanssi_out;
[i,indx_smp]=min(abs(c_CO_erot_nollasta_smp));
[i,indx_Hiu]=min(abs(c_CO_erot_nollasta_Hiu));
[i,indx_out]=min(abs(c_CO_erot_nollasta_out));
c_CO_numeerinen_smp = CO_arvauskonsentraatio(indx_smp);
c_CO_numeerinen_Hiu = CO_arvauskonsentraatio(indx_Hiu);
c_CO_numeerinen_out = CO_arvauskonsentraatio(indx_out);
disp(' ');
disp('CO: (yksittäin, suodin NB4740)');
disp('      Sample      Kiltahuone  Ulkoilma');
disp([c_CO_numeerinen_smp, c_CO_numeerinen_Hiu, c_CO_numeerinen_out]*1e6);
%disp([CO_mitattu_transmittanssi_smp,CO_mitattu_transmittanssi_Hiu,CO_mitattu_transmittanssi_out])
% Nonneh, jos toi antaa myös väärää vastausta niin meidän mittauksissa oli
% jotain ihan vinossa

h2Fig=figure (10);
set(h2Fig, 'Position', [100 100 700 600]) %rikkoo asioita
plot(CO_arvauskonsentraatio*1e6, CO_transmittanssit_exp,  "linewidth", lw,
    CO_arvauskonsentraatio*1e6, CO_transmittanssit_lin,  "linewidth", lw,
    [CO_c_min,CO_c_max]*1e6, [1,1].*CO_mitattu_transmittanssi_smp,  "linewidth", lw,
    [CO_c_min,CO_c_max]*1e6, [1,1].*CO_mitattu_transmittanssi_Hiu,  "linewidth", lw,
    [CO_c_min,CO_c_max]*1e6, [1,1].*CO_mitattu_transmittanssi_out,  "linewidth", lw)
title("Transmittanssi CO, NB4740-suotimella",'FontSize', 20)
#{
legendi=legend("Todellinen transmittanssi",
    "Transmittanssin lineaarinen approksimaatio",
    "Mitattu transmittanssi, laboratorion näyte",
    "Mitattu transmittanssi, Hiukkasen kiltahuone", 
    "Mitattu transmittanssi, ulkoilma");
set(legendi,'FontSize',18);
#}
xlabel("Pitoisuus (ppm)",'FontSize', 20);
ylabel("Transmittanssi",'FontSize', 20);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ma 5.12 muutokset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lisätty nyt veden H2O ja metaaninkin CH4 yksittäislaskut
% Kommentoitu pois minun parannukseni normalisoinnissa. koska oletan että jos 
% assaritkin ovat laskeneet arvot virheellisellä koodilla niin niin pitäisi meidänkin

% Vesi H20 ensimmäisellä suotimella NB2710
H2O_c_min = 0;
H2O_c_max = 1e-2;
H2O_arvauskonsentraatio=linspace(H2O_c_min,H2O_c_max,100);
H2O_transmittanssit_exp = sum(NB2710.*exp(-H2O.*(1e-4*p*L/(k*T))*H2O_arvauskonsentraatio).*dlambda*1e-2)/(sum(NB2710.*dlambda*1e-2));
H2O_transmittanssit_lin = 1 - ((sum(NB2710.*H2O.*H2O_arvauskonsentraatio.*dlambda*1e-4)*(1e-2*p*L/(k*T))))/(sum(NB2710.*dlambda*1e-2));
H2O_mitattu_transmittanssi_smp = (I_smp(1)/(I0_M1(1)));
H2O_mitattu_transmittanssi_Hiu = (I_Hiu(1)/(I0_M2(1)));
H2O_mitattu_transmittanssi_out = (I_out(1)/(I0_M2(1)));
c_H2O_erot_nollasta_smp = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_smp;
c_H2O_erot_nollasta_Hiu = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_Hiu;
c_H2O_erot_nollasta_out = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_out;
[i,indx_smp]=min(abs(c_H2O_erot_nollasta_smp));
[i,indx_Hiu]=min(abs(c_H2O_erot_nollasta_Hiu));
[i,indx_out]=min(abs(c_H2O_erot_nollasta_out));
c_H2O_numeerinen_smp = H2O_arvauskonsentraatio(indx_smp);
c_H2O_numeerinen_Hiu = H2O_arvauskonsentraatio(indx_Hiu);
c_H2O_numeerinen_out = H2O_arvauskonsentraatio(indx_out);
disp('H2O: (yksittäin, suodin NB2710)');
disp([c_H2O_numeerinen_smp, c_H2O_numeerinen_Hiu, c_H2O_numeerinen_out]*1e6);
h2Fig=figure (11);
set(h2Fig, 'Position', [100 100 700 600])
plot(H2O_arvauskonsentraatio*1e6, H2O_transmittanssit_exp,  "linewidth", lw,
    H2O_arvauskonsentraatio*1e6, H2O_transmittanssit_lin,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_smp,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_Hiu,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_out,  "linewidth", lw)
title("Transmittanssi H2O, NB2710-suotimella",'FontSize', 20)




% Vesi H20 viimeisellä suotimella NB5756
H2O_c_min = 0;
H2O_c_max = 1e-2;
H2O_arvauskonsentraatio=linspace(H2O_c_min,H2O_c_max,100);
H2O_transmittanssit_exp = sum(NB5756.*exp(-H2O.*(1e-4*p*L/(k*T))*H2O_arvauskonsentraatio).*dlambda*1e-2)/(sum(NB5756.*dlambda*1e-2));
H2O_transmittanssit_lin = 1 - ((sum(NB5756.*H2O.*H2O_arvauskonsentraatio.*dlambda*1e-4)*(1e-2*p*L/(k*T))))/(sum(NB5756.*dlambda*1e-2));
H2O_mitattu_transmittanssi_smp = (I_smp(7)/(I0_M1(7)));
H2O_mitattu_transmittanssi_Hiu = (I_Hiu(7)/(I0_M2(7)));
H2O_mitattu_transmittanssi_out = (I_out(7)/(I0_M2(7)));
c_H2O_erot_nollasta_smp = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_smp;
c_H2O_erot_nollasta_Hiu = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_Hiu;
c_H2O_erot_nollasta_out = H2O_transmittanssit_exp - H2O_mitattu_transmittanssi_out;
[i,indx_smp]=min(abs(c_H2O_erot_nollasta_smp));
[i,indx_Hiu]=min(abs(c_H2O_erot_nollasta_Hiu));
[i,indx_out]=min(abs(c_H2O_erot_nollasta_out));
c_H2O_numeerinen_smp = H2O_arvauskonsentraatio(indx_smp);
c_H2O_numeerinen_Hiu = H2O_arvauskonsentraatio(indx_Hiu);
c_H2O_numeerinen_out = H2O_arvauskonsentraatio(indx_out);
disp('H2O: (yksittäin, suodin NB5756)');
disp([c_H2O_numeerinen_smp, c_H2O_numeerinen_Hiu, c_H2O_numeerinen_out]*1e6);
h2Fig=figure (12);
set(h2Fig, 'Position', [100 100 700 600])
plot(H2O_arvauskonsentraatio*1e6, H2O_transmittanssit_exp,  "linewidth", lw,
    H2O_arvauskonsentraatio*1e6, H2O_transmittanssit_lin,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_smp,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_Hiu,  "linewidth", lw,
    [H2O_c_min,H2O_c_max]*1e6, [1,1].*H2O_mitattu_transmittanssi_out,  "linewidth", lw)
title("Transmittanssi H2O, NB5756-suotimella",'FontSize', 20)


% Vain veden laskeminen matriisimenetelmällä, eli ainoastaan molemmat suotimet NB2710 NB5756

I_smp_H2O=[I_smp(1);I_smp(7)];                   % Mittauksemme
I0_M1_H2O=[I0_M1(1);I0_M1(7)];
I_Hiu_H2O=[I_Hiu(1);I_Hiu(7)];
I_out_H2O=[I_out(1);I_out(7)];
I0_M2_H2O=[I0_M2(1);I0_M2(7)];
suotimet_H2O= [suotimet(:,1), suotimet(:,7)];   %suodinvektorit
absorptio_matr_H2O = absorptio_matr(:,1);        %absorptiovektorit
s_lkm_H2O=2;    %suotimia 2
% b ----------------
b_vec_smp_H2O=(1-I_smp_H2O./I0_M1_H2O) .* (sum( suotimet_H2O .*(dlambda*ones(1,s_lkm_H2O)) .*1e-2, 1))';  % Huom. transpoosi '
b_vec_Hiu_H2O=(1-I_Hiu_H2O./I0_M2_H2O) .* (sum( suotimet_H2O .*(dlambda*ones(1,s_lkm_H2O)) .*1e-2, 1))';
b_vec_out_H2O=(1-I_out_H2O./I0_M2_H2O) .* (sum( suotimet_H2O .*(dlambda*ones(1,s_lkm_H2O)) .*1e-2, 1))';
% A --------------------
sarake1_H2O = sum(absorptio_matr_H2O(:,1)*ones(1,s_lkm_H2O).*suotimet_H2O*1e-2.*(dlambda*ones(1,s_lkm_H2O)), 1)';   
A_H2O = [sarake1_H2O].*p./(k*T)*L;
c_vec_smp_H2O = A_H2O\b_vec_smp_H2O;
c_vec_Hiu_H2O = A_H2O\b_vec_Hiu_H2O;
c_vec_out_H2O = A_H2O\b_vec_out_H2O;
kaasut_ppm_H2O = [c_vec_smp_H2O, c_vec_Hiu_H2O ,c_vec_out_H2O]*1e6;
disp('H2O: Matriisimenetelmällä, käyttäen vain suotimia NB2710 NB5756');
disp(kaasut_ppm_H2O);





% Metaani CH4 suotimella NB3220
CH4_c_min = 0;
CH4_c_max = 4e-3;
CH4_arvauskonsentraatio=linspace(CH4_c_min,CH4_c_max,100);
CH4_transmittanssit_exp = sum(NB3220.*exp(-CH4.*(1e-4*p*L/(k*T))*CH4_arvauskonsentraatio).*dlambda*1e-2)/(sum(NB3220.*dlambda*1e-2));
CH4_transmittanssit_lin = 1 - ((sum(NB3220.*CH4.*CH4_arvauskonsentraatio.*dlambda*1e-4)*(1e-2*p*L/(k*T))))/(sum(NB3220.*dlambda*1e-2));
CH4_mitattu_transmittanssi_smp = (I_smp(3)/(I0_M1(3)));
CH4_mitattu_transmittanssi_Hiu = (I_Hiu(3)/(I0_M2(3)));
CH4_mitattu_transmittanssi_out = (I_out(3)/(I0_M2(3)));
c_CH4_erot_nollasta_smp = CH4_transmittanssit_exp - CH4_mitattu_transmittanssi_smp;
c_CH4_erot_nollasta_Hiu = CH4_transmittanssit_exp - CH4_mitattu_transmittanssi_Hiu;
c_CH4_erot_nollasta_out = CH4_transmittanssit_exp - CH4_mitattu_transmittanssi_out;
[i,indx_smp]=min(abs(c_CH4_erot_nollasta_smp));
[i,indx_Hiu]=min(abs(c_CH4_erot_nollasta_Hiu));
[i,indx_out]=min(abs(c_CH4_erot_nollasta_out));
c_CH4_numeerinen_smp = CH4_arvauskonsentraatio(indx_smp);
c_CH4_numeerinen_Hiu = CH4_arvauskonsentraatio(indx_Hiu);
c_CH4_numeerinen_out = CH4_arvauskonsentraatio(indx_out);
disp('CH4: (yksittäin, suodin NB3220)');
disp([c_CH4_numeerinen_smp, c_CH4_numeerinen_Hiu, c_CH4_numeerinen_out]*1e6);
h2Fig=figure (13);
set(h2Fig, 'Position', [100 100 700 600])
plot(CH4_arvauskonsentraatio*1e6, CH4_transmittanssit_exp,  "linewidth", lw,
    CH4_arvauskonsentraatio*1e6, CH4_transmittanssit_lin,  "linewidth", lw,
    [CH4_c_min,CH4_c_max]*1e6, [1,1].*CH4_mitattu_transmittanssi_smp,  "linewidth", lw,
    [CH4_c_min,CH4_c_max]*1e6, [1,1].*CH4_mitattu_transmittanssi_Hiu,  "linewidth", lw,
    [CH4_c_min,CH4_c_max]*1e6, [1,1].*CH4_mitattu_transmittanssi_out,  "linewidth", lw)
title("Transmittanssi CH4:n pitoisuuden funktiona NB3220-suotimen alueella",'FontSize', 20)
legendi=legend("Todellinen transmittanssi",
    "Transmittanssin lineaarinen approksimaatio",
    "Mitattu transmittanssi, laboratorion näyte",
    "Mitattu transmittanssi, Hiukkasen kiltahuone", 
    "Mitattu transmittanssi, ulkoilma");
set(legendi,'FontSize',18);
xlabel("Pitoisuus (ppm)",'FontSize', 20);
ylabel("Transmittanssi",'FontSize', 20);

disp('Jos pitoisuus on 0, tai määrittelyjoukon maksimi(4000,10000), on se alueen ulkopuoleella');

endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Vertailuarvoja
% H2O:  7500ppm          ulkoilma 90% suht.kosteus +5 asteessa  http://www.tis-gdv.de/tis_e/misc/klima.htm
% H2O:  6500ppm          sisäilma 30% suht.kosteus +20 asteessa --ppm = 1e6*(taulukon_arvo(g/m^3))/(18(mol/g)) / (1000l/22.4(l/mol))
% CO2:  400ppm           https://en.wikipedia.org/wiki/Carbon_dioxide
% CO:   0.2ppm           https://en.wikipedia.org/wiki/MOPITT
% CH4:  1.9ppm           https://en.wikipedia.org/wiki/Methane
% C2H4: 4ppm             http://www.inchem.org/documents/sids/sids/74851.pdf



%-----------------------------------------
% Transmissiokuvaajan piirto kaasuille ja suotimille

hFig=figure(3);
set(hFig, 'Position', [100 100 700 600]) %rikkoo asioita

q=-0.3e20;   % Täysin mielivaltainen kerroin jolla saadan mahdutettua samaan kuvaan kaksi eri asiaa
lw=2;        % linewidth

plot(lambda,H2O*q+200, "linewidth", lw,
    lambda,min((CO2), 100*(1/0.3e20))*q+200, "linewidth", lw,
    lambda,CO*q+200, "linewidth", lw,
    lambda,CH4*q+200, "linewidth", lw,
    lambda,C2H2*q+200, "linewidth", lw,
    
    [1],[1],'w',"linewidth", 0.5,
    [1],[1],'w',"linewidth", 0.5,
    [1],[1],'w',"linewidth", 0.5,
    [1],[1],'w',"linewidth", 0.5,
    
    lambda,NB2710, "linewidth", lw,
    lambda,NB3060, "linewidth", lw,
    lambda,NB3220, "linewidth", lw,
    lambda,NB3977, "linewidth", lw,
    lambda,NB4270, "linewidth", lw,
    lambda,NB4740, "linewidth", lw,
    lambda,NB5756, "linewidth", lw,
    
    [1],[1],'w',"linewidth", 0.5,
    [1],[1],'w',"linewidth", 0.5,
    
    lambda,100.*ones(1,length(lambda)),'k',"linewidth", 0.5
    )

% Vanha versio kuvaajasta
%plot(lambda,-NB2710, lambda,-NB3060,lambda,-NB3220,lambda,-NB3977,lambda,-NB4270,lambda,-NB4740,lambda,-NB5756,lambda,H2O*0.3e20, lambda,CO2*0.3e20, lambda,CO*0.3e20, lambda,CH4*0.3e20, lambda,C2H2*0.3e20)

xlabel('Aallonpituus [nm]', 'FontSize', 20)
ylabel('Transmittanssi, suotimet     Transmittanssi, kaasut', 'FontSize', 20)
title('Kaasujen ja suotimien transmittanssivertailu', 'FontSize', 20)
legendi=legend('H2O', 'CO2', 'CO', 'CH4', 'C2H2',
' ',' ',' ',' ',
'NB2710','NB3060','NB3220','NB3977','NB4270','NB4740','NB5756');
set(legendi,'FontSize',20);
set(gca,'FontSize',15)      %x-akselin tickien fonttikoko
set(gca,'ytick',[])

axis([2500,7000,0,200])


%%%% Konsentraatioiden tulostus, siirretty loppuun koska helppolukuisuus ja 
%%%% tekemiäni funktioita ei voi siirtää aiemmaksi.
disp(' ');
disp('      Sample      Kiltahuone   Ulkoilma');
disp(kaasut_ppm);
disp('Järjestys: H20, CO2, CO, CH4, C2H2');

disp(' ');
disp('CO2:');
disp('      Sample      Kiltahuone  Ulkoilma');
disp([c_CO2_numeerinen_smp, c_CO2_numeerinen_Hiu, c_CO2_numeerinen_out]*1e6);

