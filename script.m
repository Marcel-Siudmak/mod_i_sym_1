clear; clc; close all;

%% 1. Dane nominalne i parametry modelu
TzewN = -20;     % temperatura zewnętrzna [°C]
TkzN  = 30;      % temperatura powietrza nawiewanego [°C]
Tw1N  = 20;      % temperatura w pokoju 1 [°C]
Tw2N  = 15;      % temperatura w pokoju 2 [°C]
PkN   = 10000;   % nominalna moc klimatyzacji [W]
a     = 2;       % stosunek Ks1 / Ks2 oraz A1 / A2

cp   = 1000;     % ciepło właściwe powietrza [J/(kg*K)]
rop  = 1.2;      % gęstość powietrza [kg/m^3]
cpro = cp * rop; % pojemność cieplna strumienia objętościowego [J/(m^3*K)]
Ac   = 50;       % powierzchnia mieszkania [m^2]
h    = 2.5;      % wysokość pomieszczeń [m]

% Pojemności cieplne pomieszczeń
A1  = (a / (a + 1)) * Ac;
A2  = (1 / (a + 1)) * Ac;
Cv1 = cpro * A1 * h;
Cv2 = cpro * A2 * h;

% Parametry statyczne (z równań równowagi w punkcie nominalnym)
fpN = PkN / (cpro * (TkzN - Tw1N));
Ks2 = PkN / (a * (Tw1N - TzewN) + (Tw2N - TzewN));
Ks1 = a * Ks2;
K0  = (PkN - Ks1 * (Tw1N - TzewN)) / (Tw1N - Tw2N);
Kz  = (K0 * Ks2) / (K0 + Ks2);

fprintf('--- WYZNACZONE PARAMETRY MODELU ---\n');
fprintf('fpN = %.4f m^3/s | Ks1 = %.4f W/K | Ks2 = %.4f W/K | K0 = %.4f W/K\n', fpN, Ks1, Ks2, K0);
fprintf('Cv1 = %.1f J/K | Cv2 = %.1f J/K\n\n', Cv1, Cv2);

%% 2. Równania równowagi i weryfikacja punktu nominalnego
Tw1_eq = @(Tkz, Tzew, fp) (cpro * fp .* Tkz + (Ks1 + Kz) * Tzew) ./ (cpro * fp + Ks1 + Kz);
Tw2_eq = @(Tw1, Tzew)     (K0 * Tw1 + Ks2 * Tzew) ./ (K0 + Ks2);

Tw1_stat = Tw1_eq(TkzN, TzewN, fpN);
Tw2_stat = Tw2_eq(Tw1_stat, TzewN);
x0       = [Tw1N; Tw2N];

fprintf('--- WERYFIKACJA PUNKTU PRACY ---\n');
fprintf('Nominalne:  Tw1 = %.4f C, Tw2 = %.4f C\n', Tw1N, Tw2N);
fprintf('Statyczne:  Tw1 = %.4f C, Tw2 = %.4f C\n\n', Tw1_stat, Tw2_stat);

%% 3. Równania stanu dla stałego przepływu fp: dx/dt = A*x + B*u, y = C*x + D*u
% Stan: x = [Tw1; Tw2], Wejścia: u = [Tkz; Tzew], Wyjście: y = Tw2 (PV)
macierz_A = @(fp) [ -(cpro*fp + Ks1 + K0)/Cv1,   K0/Cv1;
                     K0/Cv2,                   -(K0 + Ks2)/Cv2 ];

A = macierz_A(fpN);
B = [ (cpro * fpN)/Cv1,  Ks1/Cv1;
      0,                 Ks2/Cv2 ];
C = [0, 1];
D = [0, 0];

%% 4. Punkty pracy i stałe czasowe
tab_Tzew = TzewN + [0,  10,  10,  10];
tab_Tkz  = TkzN  + [0,   0, -10, -10];
tab_Fp   = fpN   * [1,   1,   1,  0.2];
tab_Tw1  = Tw1_eq(tab_Tkz, tab_Tzew, tab_Fp);
tab_Tw2  = Tw2_eq(tab_Tw1, tab_Tzew);

etykiety = {
    '1) T_{zew}=T_{zewN}, T_{kz}=T_{kzN}, f_p=f_{pN}'
    '2) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}, f_p=f_{pN}'
    '3) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}-10, f_p=f_{pN}'
    '4) T_{zew}=T_{zewN}+10, T_{kz}=T_{kzN}-10, f_p=0.2 f_{pN}'
};

fprintf('--- PUNKTY PRACY I STAŁE CZASOWE ---\n');
for i = 1:numel(tab_Tzew)
    tau = sort(-1 ./ eig(macierz_A(tab_Fp(i))));
    fprintf('%s\n  Tw1 = %.2f C, Tw2 = %.2f C | tau1 = %.1f s, tau2 = %.1f s\n', ...
        etykiety{i}, tab_Tw1(i), tab_Tw2(i), tau(1), tau(2));
end
fprintf('\n');

%% 5. Symulacje odpowiedzi skokowych i zapis wykresów
t_skok    = 200;    % [s] chwila skoku (używana w blokach Step w Simulinku)
czas_stop = 2500;   % [s] czas symulacji (dla pkt 4 stała czasowa ~330 s)
stylLinii = {'-', '--', ':', '-.'};

% Kolejne testy: [dTkz, dTzew, dfp]
skoki = [
    0, 0, 0;
    5, 0, 0;
    0, 5, 0;
    0, 0, 0.1 * fpN
];
nazwy_testow = {
    'Test równowagi (bez skoku)'
    'Skok T_{kz} o +5 °C'
    'Skok T_{zew} o +5 °C'
    'Skok f_p o +10% f_{pN}'
};

tytuly_osi = {'T_{w1}', 'T_{w2}', '\DeltaT_{w1} = T_{w1} - T_{w1}(0)', '\DeltaT_{w2} = T_{w2} - T_{w2}(0)'};
opisy_y    = {'Temperatura [°C]', 'Temperatura [°C]', '\DeltaT [°C]', '\DeltaT [°C]'};

if ~exist('./plot', 'dir')
    mkdir('./plot');
end

for j = 1:size(skoki, 1)
    dTkz  = skoki(j, 1);
    dTzew = skoki(j, 2);
    dfp   = skoki(j, 3);

    fig = figure('Name', nazwy_testow{j});
    tl  = tiledlayout(2, 2);
    if j == 1
        title(tl, nazwy_testow{j});
    else
        title(tl, sprintf('%s w t = %g s', nazwy_testow{j}, t_skok));
    end

    ax = gobjects(1, 4);
    for k = 1:4
        ax(k) = nexttile;
        hold(ax(k), 'on');
        grid(ax(k), 'on');
        title(ax(k), tytuly_osi{k});
        xlabel(ax(k), 'Czas [s]');
        ylabel(ax(k), opisy_y{k});
    end

    for i = 1:numel(tab_Tzew)
        % Zmienne przekazywane do modelu Simulink ('schemat')
        Tkz_0  = tab_Tkz(i);
        Tzew_0 = tab_Tzew(i);
        fp_0   = tab_Fp(i);
        Tw1_0  = tab_Tw1(i);
        Tw2_0  = tab_Tw2(i);

        out = sim('schemat', 'StopTime', num2str(czas_stop), 'MaxStep', '1');
        t   = out.simout.Time;
        Tw1 = out.simout.Data(:, 1);
        Tw2 = out.simout.Data(:, 2);

        sygnaly = {Tw1, Tw2, Tw1 - Tw1(1), Tw2 - Tw2(1)};
        for k = 1:4
            plot(ax(k), t, sygnaly{k}, 'LineWidth', 1.5, 'LineStyle', stylLinii{i});
        end

        % Porównanie wartości ustalonej z symulacji z rozwiązaniem statycznym
        Tw1_k = Tw1_eq(Tkz_0 + dTkz, Tzew_0 + dTzew, fp_0 + dfp);
        Tw2_k = Tw2_eq(Tw1_k, Tzew_0 + dTzew);
        fprintf('%s | pkt %d: sym [%.3f, %.3f] | obl [%.3f, %.3f] | dT [%.3f, %.3f]\n', ...
            nazwy_testow{j}, i, Tw1(end), Tw2(end), Tw1_k, Tw2_k, ...
            Tw1(end) - Tw1(1), Tw2(end) - Tw2(1));
    end

    lgd = legend(ax(1), etykiety, 'NumColumns', 2);
    lgd.Layout.Tile = 'south';

    exportgraphics(fig, sprintf('./plot/wykres_%d.png', j), 'Resolution', 300);
end
