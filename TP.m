%% ECHANTILLON Yoan & BLANC Nicolas
clear;
close all;
clc;
%% Initialisation des paramètres
fech=10E3;
Ds=1E3;
Fse=fech/Ds;
Ts = 1e-3;

M=4;
n_b=log2(M);
Ak=[(-1-1j)/sqrt(2); (-1+1j)/sqrt(2); (1-1j)/sqrt(2); (1+1j)/sqrt(2)];

Ns=5000;
Nc=Ns*Fse;
Nfft=512;

%freq porteuse
f0 = 2.5E3;

G=rcosdesign(0.35,8,4,'sqrt');  %filtre de mise en forme en racin de cosinus surélevé

Eg = 0; % Energie du filtre de mise en forme -> somme des modules au carré 
for i=1:length(G)
    Eg = Eg + G(i)^2;
end

sigA2 = 1; % Variance théorique des symboles -> calcul a partir de la formule avec E(X)²

eb_n0_dB = 0:0.5:10; % Liste des Eb/N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit réel sur freq porteuse

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB (résultats)
Pb = qfunc ( sqrt (2* eb_n0 ) ) ; % Tableau des probabilités d’erreurs théoriques = 0.5*erfc(sqrt(eb_n0))

for j = 1: length(eb_n0)
    bit_error = 0;
    bit_count = 0;
    while bit_error < 100
%% Emetteur %%

    Sb=randi([0,3],1,Ns); %génère Ns échantillons aléatoires entre 0 et 3 (00,01,10,11)
    %1 échantillon = 2 bits           
                      
    Ss = pskmod(Sb,M,pi/4,'gray'); %bit->symbole

    Ss2=upsample(Ss,Fse); %suréchantillonnage
    
    Sl=conv2(G,Ss2);   
    
    t_port=(0:length(Sl)-1)/fech;
    port_cplx = exp(1i*2*pi*f0*t_port);

    S = real(Sl .*port_cplx); % enveloppe cplx du signal
%% CANAL
    alpha0 = 1;
    nl = sqrt(sigma2(j)/4)* (randn(size(S)));  %bruit blanc réel
    y = alpha0 * S + nl;
    
%% ETAGE RF->BDB

    freq_error = 0;

    port = 2*cos(2*pi*(f0+freq_error)*t_port);
    port2 = 2*cos(2*pi*(f0+freq_error)*t_port + pi/2);
    
    yI = y .*port;
    yQ = y .*port2;

    yl = yI + 1i*yQ;
    
%% RECEPTEUR

    Ga = G; %filtre adapté on peut aussi écrire fliplr(G) mais comme c'est symétrique c'est identique
    Rg = conv2(G,Ga); %Autocorrélation entre le filtre G et le filtre adapaté Ga
    
    retard = 0;
    max = Rg(1);
    for i=2:length(Rg)  %calcul du retard lié aux filtres
        if (Rg(i) > max)
            retard = i;
            max = Rg(i);
        end
    end
    
    rl = conv2(Ga, yl);

    rln = rl(retard:Fse:length(rl)); %sous-echantillonnage
                                     %*exp(1i*pi/6) pour tout décaler de pi/6 
 
    bn = pskdemod(rln,M,pi/4,'gray'); % décision & Symbole -> bit
    
    Sb2 = zeros(1, n_b*length(Sb)); % 0,1,2,3 -> 00,01,10,11 pour évaluer bit à bit
    for i=1:1:length(Sb)
        if (Sb(i) == 0)
            Sb2(2*i-1) = 0; Sb2(2*i) = 0;
        elseif (Sb(i) == 1)
            Sb2(2*i-1) = 0; Sb2(2*i) = 1;
        elseif (Sb(i) == 2)
            Sb2(2*i-1) = 1; Sb2(2*i) = 0;
        else
            Sb2(2*i-1) = 1; Sb2(2*i) = 1;
        end
    end
    
    bn2 = zeros(1, n_b*length(bn)); % 0,1,2,3 -> 00,01,10,11 pour évaluer bit à bit
    for i=1:1:length(bn)
        if (bn(i) == 0)
            bn2(2*i-1) = 0; bn2(2*i) = 0;
        elseif (bn(i) == 1)
            bn2(2*i-1) = 0; bn2(2*i) = 1;
        elseif (bn(i) == 2)
            bn2(2*i-1) = 1; bn2(2*i) = 0;
        else
            bn2(2*i-1) = 1; bn2(2*i) = 1;
        end
    end

    for i =1:Ns*n_b
        if(Sb2(i) ~= bn2(i))
            bit_error = bit_error + 1;
        end
        bit_count = bit_count + 1;
    end
    end

    TEB(j) = bit_error/bit_count;
end


%% Affichage des résultats
Te = 1/fech;
% t_plot = linspace(0,20*Ts-Te, 200);
% 
% 
% figure(1);
% plot(t_plot,real(Sl(1:200)),'b');
% 
% hold on;
% plot(t_plot,real(rl(1:200)),'r');
% xlabel("temps en s");
% ylabel("partie réelle de s_l(t) et r_l(t)");
% title("s_l(t) en bleu et r_l(t) en rouge");

scatterplot(rln);

figure(4);
semilogy(eb_n0_dB,TEB,'b');
hold on
semilogy(eb_n0_dB,Pb,'r');
xlabel("E_b/N_0 en dB");
ylabel("log(TEB)");
title("évolution du TEB en fonction du SNR");



freq = linspace(-fech/2,fech/2, Nfft);

[DSP_exp,f]=pwelch(S,ones(1,Nfft),0,Nfft,fech,'centered');
Tf_S = fft(S,Nfft);
DSP_th = abs((Tf_S)).^2;

figure(2);
semilogy(freq,fftshift(DSP_exp),'b'); %facteur 5000 d'écart
hold on;

semilogy(freq,fftshift(real(DSP_th)),'r');
xlabel("axe fréquentiel (Hz)");
ylabel("DSP (s(t))");
title("DSP Welch en bleu DSP th en rouge");
