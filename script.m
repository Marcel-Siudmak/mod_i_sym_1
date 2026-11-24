clear; clc; close all;

%% 1. Dane nominalne i wyznaczenie parametrów modelu
% Wartości nominalne temperatur [st. C] i mocy [W]


TzewN = -20;
TkzN  = 30;
Tw1N  = 20;
Tw2N  = 15;
PkN   = 10000;     % 10 kW
a     = 2;        % stosunek Ks1 / Ks2

% Stałe fizyczne powietrza i wymiary mieszkania
cp   = 1000;    % ciepło właściwe powietrza [J/(kg*K)]
rop  = 1.2;     % gęstość powietrza [kg/m^3]
Ac   = 50;      % całkowita powierzchnia mieszkania [m^2]
h    = 2.5;     % wysokość pomieszczeń [m]

%

% Podzial powierzchni (proporcjonalnie do a = 2) i pojemnosci cieplne
A1  = (a / (a + 1)) * Ac;
A2  = (1 / (a + 1)) * Ac;
Cv1 = cp * rop * A1 * h;
Cv2 = cp * rop * A2 * h;

% Obliczenie parametrow statycznych
fpN = PkN / (cp * rop * (TkzN - Tw1N));
Ks2 = PkN / (a * (Tw1N - TzewN) + (Tw2N - TzewN));
Ks1 = a * Ks2;
K0  = (PkN - Ks1 * (Tw1N - TzewN)) / (Tw1N - Tw2N);



fprintf('--- WYZNACZONE PARAMETRY MODELU ---\n');
fprintf('fpN = %.4f m^3/s\n', fpN);
fprintf('Ks1 = %.4f W/K\n', Ks1);
fprintf('Ks2 = %.4f W/K\n', Ks2);
fprintf('K0  = %.4f W/K\n', K0);
fprintf('Cv1 = %.1f J/K, Cv2 = %.1f J/K\n\n', Cv1, Cv2);

%% 2. Obliczenie punktu pracy – rozwiązanie układu równań statycznych
fp   = fpN;
Tkz  = TkzN;
Tzew = TzewN;
cpro = cp*rop;   % moc cieplna na 1 K różnicy temperatur [W/K]

Kz = (K0 * Ks2) / (K0 + Ks2);
Tw1_stat = (cp * rop * fp * Tkz + (Ks1 + Kz) * Tzew) / (cp * rop * fp + Ks1 + Kz);
Tw2_stat = (K0 * Tw1_stat + Ks2 * Tzew) / (K0 + Ks2);


%% 3. Wartości początkowe (warunki początkowe dla symulacji / integratorów)
Tw1_0 = Tw1N;                 % wartość początkowa stanu Tw1 w chwili t=0 [st. C]
Tw2_0 = Tw2N;                 % wartość początkowa stanu Tw2 w chwili t=0 [st. C]
x0    = [Tw1_0; Tw2_0];       % wektor warunków początkowych stanu x(0)


%% 4. Zastosowanie równań stanu: dx/dt = A*x + B*u, y = C*x + D*u
% Stan: x = [Tw1; Tw2], Wejścia: u = [Tkz; Tzew], Wyjście: y = Tw2 (PV)
A = [ -(cp*rop*fp + Ks1 + K0)/Cv1,   K0/Cv1;
    K0/Cv2,                    -(K0 + Ks2)/Cv2 ];

B = [ (cp*rop*fp)/Cv1,  Ks1/Cv1;
    0,               Ks2/Cv2 ];

C = [0, 1];                   % wyjście regulowane PV = Tw2
D = [0, 0];


%% 5. Sprawdzenie poprawności obliczeń (punkt pracy == wartości nominalne)
fprintf('--- WERYFIKACJA PUNKTU PRACY ---\n');
fprintf('Wartości nominalne:   Tw1N = %.4f C, Tw2N = %.4f C\n', Tw1N, Tw2N);
fprintf('Równania statyczne:   Tw1  = %.4f C, Tw2  = %.4f C\n', Tw1_stat, Tw2_stat);
fprintf('Warunki początkowe:   x1(0)= %.4f C, x2(0)= %.4f C\n', x0(1), x0(2));


%% 6. Zadanie 3: odpowiedzi skokowe w różnych punktach pracy
t_skok    = 200;      % [s] chwila skoku (taki sam Step time w blokach Step)
czas_stop = 2500;     % [s] w punkcie 4 stała czasowa to ok. 330 s
stylLinii = {'-','--',':','-.'};

% --- równowaga dla zadanych wejść
Kz   = K0*Ks2/(K0+Ks2);
Tw1f = @(Tkz,Tzew,fp) (cp*rop*fp.*Tkz + (Ks1+Kz)*Tzew)./(cp*rop*fp + Ks1 + Kz);
Tw2f = @(Tw1,Tzew)    (K0*Tw1 + Ks2*Tzew)./(K0+Ks2);

% --- tabela punktów pracy
tab_Tzew = TzewN + [0   10   10   10];
tab_Tkz  = TkzN  + [0    0  -10  -10];
tab_Fp   = fpN   * [1    1    1   0.2];
tab_Tw1  = Tw1f(tab_Tkz, tab_Tzew, tab_Fp);
tab_Tw2  = Tw2f(tab_Tw1, tab_Tzew);
n = numel(tab_Tzew);

etyk = {'1) T_{zew}=T_{zewN}, T_{kz}=T_{kzN}, f_p=f_{pN}', ...
    '2) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}, f_p=f_{pN}', ...
    '3) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}-10, f_p=f_{pN}', ...
    '4) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}-10, f_p=0.2 f_{pN}'};

% --- stałe czasowe w każdym punkcie
fprintf('\n--- PUNKTY PRACY I STAŁE CZASOWE ---\n');
for i = 1:n
    Ai  = [-(cp*rop*tab_Fp(i)+Ks1+K0)/Cv1, K0/Cv1; K0/Cv2, -(K0+Ks2)/Cv2];
    tau = sort(-1./eig(Ai));
    fprintf('%s\n  Tw1=%.2f, Tw2=%.2f | tau1=%.1f s, tau2=%.1f s\n', ...
        etyk{i}, tab_Tw1(i), tab_Tw2(i), tau(1), tau(2));
end

% --- skoki: [dTkz dTzew dfp], ten sam skok w każdym punkcie pracy
skoki = [0 0 0;  5 0 0;  0 5 0;  0 0 0.1*fpN];
nazwy = {'Test równowagi (bez skoku)', ...
    'Skok T_{kz} o +5 °C', ...
    'Skok T_{zew} o +5 °C', ...
    'Skok f_p o +10% f_{pN}'};

if ~exist('./plot', 'dir')
    mkdir('./plot');
end


for j = 1:size(skoki,1)
    dTkz = skoki(j,1);  dTzew = skoki(j,2);  dfp = skoki(j,3);
    figure('Name', nazwy{j});
    tl = tiledlayout(2,2);
    if j == 1
        title(tl, nazwy{j});
    else
        title(tl, [nazwy{j} ' w t = ' num2str(t_skok) ' s']);
    end
    ax = [nexttile nexttile nexttile nexttile]; hold(ax,'on');

    for i = 1:n
        % wartości dla Simulinka 
        Tkz_0 = tab_Tkz(i);  Tzew_0 = tab_Tzew(i);  fp_0 = tab_Fp(i);
        Tw1_0 = tab_Tw1(i);  Tw2_0  = tab_Tw2(i);

        out = sim('schemat', 'StopTime', num2str(czas_stop), 'MaxStep', '1');
        t = out.simout.Time;  y = out.simout.Data;

        plot(ax(1), t, y(:,1), 'LineWidth', 1.5, 'LineStyle', stylLinii{i});
        plot(ax(2), t, y(:,2), 'LineWidth', 1.5, 'LineStyle', stylLinii{i});
        plot(ax(3), t, y(:,1)-y(1,1), 'LineWidth', 1.5, 'LineStyle', stylLinii{i});
        plot(ax(4), t, y(:,2)-y(1,2), 'LineWidth', 1.5, 'LineStyle', stylLinii{i});

        % kontrola: wartość ustalona z symulacji vs z równań statycznych
        Tw1_k = Tw1f(Tkz_0+dTkz, Tzew_0+dTzew, fp_0+dfp);
        Tw2_k = Tw2f(Tw1_k, Tzew_0+dTzew);
        fprintf('%s | punkt %d: sym Tw1=%.3f Tw2=%.3f | obl Tw1=%.3f Tw2=%.3f | dTw1=%.3f dTw2=%.3f\n', ...
            nazwy{j}, i, y(end,1), y(end,2), Tw1_k, Tw2_k, ...
            y(end,1)-y(1,1), y(end,2)-y(1,2));
    end

    tyt = {'T_{w1}','T_{w2}','\DeltaT_{w1} = T_{w1} - T_{w1}(0)','\DeltaT_{w2} = T_{w2} - T_{w2}(0)'};
    for k = 1:4
        grid(ax(k),'on'); xlabel(ax(k),'Czas [s]'); title(ax(k),tyt{k});
        if k <= 2, ylabel(ax(k),'Temperatura [°C]'); else, ylabel(ax(k),'\DeltaT [°C]'); end
    end
    lgd = legend(ax(1), etyk, 'NumColumns', 2);
    lgd.Layout.Tile = 'south';

    exportgraphics(gcf, sprintf('./plot/wykres_%d.png', j), 'Resolution', 300);
end