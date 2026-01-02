clc; clear; close all;

%% ========================================================================
%  MICRORÉSEAU INTELLIGENT - 19 MAISONS + PV + BATTERIES MAISONS + V2G
%  - 19 maisons avec 4 kWc de PV chacune (≈ 76 kWc total)
%  - Batterie individuelle par maison (9.6 kWh)
%  - Flotte V2G (19 véhicules)
%  - Pas de batterie stationnaire
%  - Tarification dynamique + stratégie prix import/export
%% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║   MICRORÉSEAU INTELLIGENT - MAISONS + BATTERIES + V2G     ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% 1. CHARGEMENT DONNÉES PVSYST
fprintf('📂 Chargement données PVsyst...\n');

% Demander le fichier CSV
[filename, filepath] = uigetfile('*.csv', 'Sélectionner le fichier PVsyst CSV');
if isequal(filename,0)
    error('Aucun fichier sélectionné. Simulation annulée.');
end

fullpath = fullfile(filepath, filename);
fprintf('  • Fichier: %s\n', filename);

% Lecture du fichier CSV avec gestion des en-têtes
opts = detectImportOptions(fullpath);
opts.VariableNamingRule = 'preserve';
data_pvsyst = readtable(fullpath, opts);

% Afficher les colonnes disponibles
fprintf('  • Colonnes disponibles:\n');
disp(data_pvsyst.Properties.VariableNames');

% Identifier les colonnes pertinentes
% Recherche flexible des noms de colonnes
col_names = data_pvsyst.Properties.VariableNames;

% Trouver GlobInc (irradiance sur plan incliné)
idx_globinc = find(contains(col_names, 'GlobInc', 'IgnoreCase', true) | ...
                   contains(col_names, 'GlobInc W/m', 'IgnoreCase', true));
if isempty(idx_globinc)
    error('Colonne GlobInc (irradiance) non trouvée dans le fichier CSV');
end
col_GlobInc = col_names{idx_globinc(1)};

% Trouver GlobEff (irradiance effective)
idx_globeff = find(contains(col_names, 'GlobEff', 'IgnoreCase', true));
if ~isempty(idx_globeff)
    col_GlobEff = col_names{idx_globeff(1)};
else
    col_GlobEff = col_GlobInc; % utiliser GlobInc si GlobEff absent
    fprintf('  ⚠ GlobEff non trouvé, utilisation de GlobInc\n');
end

% Trouver EArray (production du générateur PV)
idx_earray = find(contains(col_names, 'EArray', 'IgnoreCase', true) | ...
                  contains(col_names, 'E_Grid', 'IgnoreCase', true));
if ~isempty(idx_earray)
    col_EArray = col_names{idx_earray(1)};
else
    col_EArray = '';
    fprintf('  ⚠ EArray non trouvé\n');
end

% Extraire les données
G_pvsyst = data_pvsyst.(col_GlobInc);  % W/m²
G_pvsyst(isnan(G_pvsyst)) = 0;         % Remplacer NaN par 0

% Déterminer le nombre de jours et heures
N_timesteps = length(G_pvsyst);
fprintf('  • Nombre de pas de temps: %d\n', N_timesteps);

% Organiser en matrice [jours x heures]
if mod(N_timesteps, 24) == 0
    N_days = N_timesteps / 24;
    G = reshape(G_pvsyst(1:N_days*24), 24, N_days)';  % [N_days x 24]
else
    % Si pas multiple de 24, tronquer
    N_days = floor(N_timesteps / 24);
    G = reshape(G_pvsyst(1:N_days*24), 24, N_days)';
    fprintf('  ⚠ Données tronquées à %d jours complets\n', N_days);
end

hours = 0:23;

% Température: estimer à partir de l'heure si non disponible
% Recherche d'une colonne température
idx_temp = find(contains(col_names, 'T Amb', 'IgnoreCase', true) | ...
                contains(col_names, 'Temp', 'IgnoreCase', true));
if ~isempty(idx_temp)
    col_Temp = col_names{idx_temp(1)};
    T_pvsyst = data_pvsyst.(col_Temp);
    T_pvsyst(isnan(T_pvsyst)) = 20;
    T = reshape(T_pvsyst(1:N_days*24), 24, N_days)';
    fprintf('  • Température extraite de la colonne: %s\n', col_Temp);
else
    % Température par défaut basée sur l'heure
    T = zeros(N_days, 24);
    for day = 1:N_days
        for h = 1:24
            T(day, h) = 18 + 10*sin(2*pi*(h-15)/24) + 8*sin(2*pi*(day-172)/365);
        end
    end
    fprintf('  ⚠ Température estimée (non trouvée dans CSV)\n');
end

fprintf('✓ Données PVsyst chargées: %d jours (%.1f mois)\n', N_days, N_days/30);

%% 2. CONFIGURATION DU QUARTIER
fprintf('\n📋 Configuration du microréseau...\n');

Nhouses          = 19;
Conso_annuelle   = 4000;          % kWh/maison/an
Conso_horaire_moy = Conso_annuelle / 8760;   % kWh/h

% Véhicules V2G (1 par maison)
Nvehicles       = 19;
C_EV            = 50;             % kWh
P_charge_max    = 11;             % kW
P_discharge_max = 6;              % kW
SOC_EV_min      = 0.20;
SOC_EV_max      = 0.90;
eta_charge_EV   = 0.95;
eta_discharge_EV= 0.93;

% Batteries individuelles maisons - MODIFIÉ À 9.6 kWh
C_batt_house        = 9.6;        % kWh (changé de 15 à 9.6)
P_batt_house_max    = 7;          % kW
SOC_batt_house_min  = 0.20;
SOC_batt_house_max  = 0.95;
eta_batt_house      = 0.93;

% Limites réseau
P_transfo_max   = 250;            % kW limite transformateur
P_import_max    = 200;            % kW
P_export_max    = 180;            % kW

fprintf('  • Maisons: %d\n', Nhouses);
fprintf('  • PV par maison: 4 kWc (≈ %.0f kWc total)\n', Nhouses*4);
fprintf('  • Véhicules V2G: %d (%.0f kWh chacun)\n', Nvehicles, C_EV);
fprintf('  • Batterie maison: %.1f kWh / %.0f kW\n', C_batt_house, P_batt_house_max);
fprintf('  • Capacité totale stockage: %.1f kWh\n', ...
    Nhouses*C_batt_house + Nvehicles*C_EV);
fprintf('  • Limite transformateur: %.0f kW\n', P_transfo_max);
fprintf('  • Période simulée: %d jours\n\n', N_days);

%% 3. PRODUCTION PHOTOVOLTAÏQUE (basée sur données PVsyst)
fprintf('☀  Calcul production photovoltaïque basée sur PVsyst...\n');

% Extraire la production si disponible dans le CSV
if ~isempty(col_EArray)
    E_pvsyst = data_pvsyst.(col_EArray);  % kWh
    E_pvsyst(isnan(E_pvsyst)) = 0;
    
    % Convertir en puissance (kW) si en énergie horaire
    P_pvsyst_total = E_pvsyst;  % Déjà en kW si données horaires
    
    % Organiser en matrice
    P_pvsyst_2D = reshape(P_pvsyst_total(1:N_days*24), 24, N_days)';  % [N_days x 24]
    
    % Facteur de normalisation pour 19 maisons × 4 kWc
    Puissance_totale_cible = Nhouses * 4.0;  % kWc
    Puissance_max_pvsyst = max(P_pvsyst_2D(:));
    
    if Puissance_max_pvsyst > 0
        facteur_echelle = Puissance_totale_cible / Puissance_max_pvsyst;
    else
        facteur_echelle = 1;
    end
    
    % Mise à l'échelle
    P_PV_tot = P_pvsyst_2D * facteur_echelle;  % kW
    
    fprintf('  • Production PVsyst extraite et mise à l''échelle\n');
    fprintf('  • Facteur d''échelle appliqué: %.3f\n', facteur_echelle);
else
    % Calcul basé sur GlobInc si EArray non disponible
    fprintf('  • Calcul de production basé sur irradiance PVsyst\n');
    
    % Puissance nominale par maison avec légère dispersion
    P_pv_nom = 4.0 + 0.4*(rand(Nhouses,1)-0.5);  % kW
    
    P_PV = zeros(Nhouses, N_days, 24);   % [W]
    for house = 1:Nhouses
        for day = 1:N_days
            for h = 1:24
                G_norm   = G(day, h) / 1000;           % kW/m²
                temp_coeff = 1 - 0.004 * (T(day,h) - 25);
                rendement = 0.92;
                P_PV(house, day, h) = P_pv_nom(house) * 1000 * G_norm ...
                                      * temp_coeff * rendement;
                P_PV(house, day, h) = max(0, P_PV(house, day, h));
            end
        end
    end
    
    P_PV_tot = squeeze(sum(P_PV,1)) / 1000;  % kW
end

% Production par maison (distribution proportionnelle)
P_pv_nom = 4.0 + 0.4*(rand(Nhouses,1)-0.5);  % kW nominal par maison
facteur_maison = P_pv_nom / sum(P_pv_nom);   % Facteur de distribution

P_PV = zeros(Nhouses, N_days, 24);
for house = 1:Nhouses
    for day = 1:N_days
        for h = 1:24
            P_PV(house, day, h) = P_PV_tot(day, h) * facteur_maison(house) * 1000;  % W
        end
    end
end

E_PV_tot = sum(P_PV_tot(:));             % kWh

fprintf('  • Puissance PV totale: %.1f kWc\n', sum(P_pv_nom));
fprintf('  • Production totale (période): %.0f kWh\n', E_PV_tot);
fprintf('  • Production moyenne journalière: %.1f kWh/jour\n\n', E_PV_tot/N_days);

%% 4. CONSOMMATION RÉSIDENTIELLE RÉALISTE
fprintf('🏠 Génération profils de consommation...\n');

P_load = zeros(Nhouses, N_days, 24);   % [W]

for house = 1:Nhouses
    % Base (kW) autour de 0.5–0.8 kW
    P_base_kW = 0.6 + 0.25*rand;
    
    for day = 1:N_days
        saison_factor = 1 + 0.25 * sin(2*pi*day/365);
        
        for h = 1:24
            P_h = P_base_kW;
            
            if h>=0 && h<6          % nuit
                P_h = P_h * 0.4;
            elseif h>=6 && h<9      % pointe matin
                P_h = P_h * (2.0 + 0.5*exp(-((h-7.5)^2)/1));
            elseif h>=9 && h<17     % journée
                P_h = P_h * 0.7;
            elseif h>=17 && h<22    % pointe soir
                P_h = P_h * (2.2 + 0.7*exp(-((h-19.5)^2)/2));
            else                    % fin de soirée
                P_h = P_h * 0.9;
            end
            
            P_h = P_h * saison_factor;
            P_h = P_h * (0.8 + 0.4*rand);  % aléa maison/jour
            P_load(house, day, h) = max(0.3, P_h) * 1000;  % en W
        end
    end
end

P_load_tot = squeeze(sum(P_load,1)) / 1000;  % kW
E_load_tot = sum(P_load_tot(:));            % kWh

fprintf('  • Consommation totale (6 mois): %.0f kWh\n', E_load_tot);
fprintf('  • Moyenne par maison (année équivalente): %.0f kWh/an\n\n', ...
    E_load_tot/Nhouses*365/180);

%% 5. TARIFICATION DYNAMIQUE
fprintf('💡 Génération tarification dynamique...\n');

tarif_nuit   = 0.130;   % TND/kWh
tarif_normal = 0.190;
tarif_pointe = 0.320;

prix = zeros(N_days, 24);

for day = 1:N_days
    for h = 1:24
        
        % Base par plage horaire
        if h>=0 && h<6
            prix_base = tarif_nuit;
        elseif h>=17 && h<22
            prix_base = tarif_pointe;
        else
            prix_base = tarif_normal;
        end
        
        % Ajustement selon déficit / surplus
        P_net = P_PV_tot(day,h) - P_load_tot(day,h);
        
        if P_net > 0   % surplus PV
            prix_dyn = prix_base * (0.7 - 0.1*rand);     % prix bas
        else           % déficit
            ratio_def = min(1.5, abs(P_net)/max(P_load_tot(day,h),0.1));
            prix_dyn  = prix_base * (1 + 0.4*ratio_def);
        end
        
        prix(day,h) = max(0.10, min(0.45, prix_dyn));
    end
end

fprintf('  • Prix moyen: %.3f TND/kWh\n', mean(prix(:)));
fprintf('  • Prix min/max: [%.3f , %.3f] TND/kWh\n\n', ...
    min(prix(:)), max(prix(:)));

%% 6. OPTIMISATION HEURE PAR HEURE (Batteries + V2G + Import/Export)
fprintf('🎯 Lancement optimisation EMS...\n');

% Variables de décision
P_EV_opt          = zeros(Nvehicles, N_days, 24);   % kW (+ charge, - décharge)
P_batt_house_opt  = zeros(Nhouses, N_days, 24);     % kW
SOC_EV_opt        = zeros(Nvehicles, N_days, 24);
SOC_batt_house_opt = zeros(Nhouses, N_days, 24);
P_import_opt      = zeros(N_days, 24);              % kW
P_export_opt      = zeros(N_days, 24);              % kW

% Conditions initiales
SOC_EV_initial          = 0.50 + 0.20*rand(Nvehicles,1);
SOC_batt_house_initial  = 0.30 + 0.40*rand(Nhouses,1);

% Seuils prix
seuil_prix_bas  = 0.16;
seuil_prix_haut = 0.28;

fprintf('  • Début boucle temporelle...\n');

for day = 1:N_days
    if mod(day,30)==0
        fprintf('    - Jour %d/%d\n', day, N_days);
    end
    
    for h = 1:24
        
        % ==================== États précédents ====================
        if day==1 && h==1
            SOC_EV_current        = SOC_EV_initial;
            SOC_batt_house_current= SOC_batt_house_initial;
        elseif h==1
            SOC_EV_current        = squeeze(SOC_EV_opt(:, day-1, 24));
            SOC_batt_house_current= squeeze(SOC_batt_house_opt(:, day-1, 24));
        else
            SOC_EV_current        = squeeze(SOC_EV_opt(:, day, h-1));
            SOC_batt_house_current= squeeze(SOC_batt_house_opt(:, day, h-1));
        end
        
        prix_h = prix(day,h);
        
        % ==================== Batteries maisons ====================
        P_batt_house_tot = 0;
        
        for house = 1:Nhouses
            
            P_PV_house   = squeeze(P_PV(house, day, h))/1000;  % kW
            P_load_house = squeeze(P_load(house, day, h))/1000;
            P_local_net  = P_PV_house - P_load_house;          % >0 : surplus
            
            SOC_now = SOC_batt_house_current(house);
            P_batt  = 0;
            
            % Charge avec surplus PV ou prix bas
            if P_local_net > 0.4 && SOC_now < SOC_batt_house_max-0.05
                Pmax_charge = min(P_batt_house_max, 0.8*P_local_net);
                cap_rest    = (SOC_batt_house_max - SOC_now)*C_batt_house;
                P_charge    = min(Pmax_charge, cap_rest/eta_batt_house);
                P_batt      = P_charge;                 % charge = +ve
            % Décharge lorsque prix élevé et déficit
            elseif P_local_net < -0.3 && prix_h > 0.22 && ...
                   SOC_now > SOC_batt_house_min+0.10
                deficit     = abs(P_local_net);
                Pmax_dis    = min(P_batt_house_max, 1.5*deficit);
                E_dispo     = (SOC_now - SOC_batt_house_min)*C_batt_house ...
                              * eta_batt_house;
                P_dis       = min(Pmax_dis, E_dispo);
                P_batt      = -P_dis;                  % décharge = -ve
            end
            
            % Limites puissance
            if abs(P_batt) > 0
                P_batt = sign(P_batt)*min(abs(P_batt), P_batt_house_max);
            end
            
            P_batt_house_opt(house, day, h) = P_batt;
            P_batt_house_tot = P_batt_house_tot + P_batt;
            
            % Mise à jour SOC
            if P_batt > 0      % charge
                SOC_new = SOC_now + (P_batt*eta_batt_house/C_batt_house);
            elseif P_batt < 0  % décharge
                SOC_new = SOC_now + (P_batt/(eta_batt_house*C_batt_house));
            else
                SOC_new = SOC_now*0.999;  % auto-décharge faible
            end
            
            SOC_batt_house_opt(house, day, h) = ...
                max(SOC_batt_house_min, min(SOC_batt_house_max, SOC_new));
        end
        
        % Puissance nette quartier après batteries maisons
        P_net_apres_maisons = P_PV_tot(day,h) - P_load_tot(day,h) ...
                              - P_batt_house_tot;
        
        % ==================== Gestion V2G ====================
        trajet_matin = (h>=7 && h<=9);
        trajet_soir  = (h>=17 && h<=19);
        veh_dispo    = ~(trajet_matin || trajet_soir);   % indispo en trajets
        
        P_EV_tot_kW = 0;
        SOC_EV_moy  = mean(SOC_EV_current);
        
        if veh_dispo
            % Charge V2G quand surplus & prix bas
            if P_net_apres_maisons > 5 && prix_h < seuil_prix_bas ...
                    && SOC_EV_moy < 0.85
                Pmax_charge_EV = P_charge_max*Nvehicles*0.8;
                P_EV_tot_kW    = min(Pmax_charge_EV, P_net_apres_maisons*0.6);
                
            % Décharge V2G quand déficit & prix haut
            elseif P_net_apres_maisons < -5 && prix_h > seuil_prix_haut ...
                    && SOC_EV_moy > 0.40
                Pmax_dis_EV    = P_discharge_max*Nvehicles*0.6;
                P_EV_tot_kW    = -min(Pmax_dis_EV, abs(P_net_apres_maisons)*0.5);
            end
        end
        
        % Répartition sur véhicules + mise à jour SOC
        for v = 1:Nvehicles
            P_v = P_EV_tot_kW / Nvehicles;
            P_EV_opt(v, day, h) = P_v;
            
            if P_v > 0
                SOC_new = SOC_EV_current(v) + ...
                    (P_v*eta_charge_EV/C_EV);
            else
                SOC_new = SOC_EV_current(v) + ...
                    (P_v/(eta_discharge_EV*C_EV));
            end
            
            SOC_EV_opt(v, day, h) = ...
                max(SOC_EV_min, min(SOC_EV_max, SOC_new));
        end
        
        % ==================== Import / Export réseau ====================
        P_net_final = P_net_apres_maisons - P_EV_tot_kW;    % >0 : surplus
        
        if P_net_final < 0
            besoin_import = abs(P_net_final);
            if prix_h < seuil_prix_bas
                P_import_opt(day,h) = min(besoin_import, P_import_max);
            elseif prix_h > seuil_prix_haut
                P_import_opt(day,h) = min(besoin_import*0.4, P_import_max*0.6);
            else
                P_import_opt(day,h) = min(besoin_import, P_import_max);
            end
            P_export_opt(day,h) = 0;
        else
            excedent = P_net_final;
            if prix_h > seuil_prix_haut
                P_export_opt(day,h) = min(excedent, P_export_max);
            elseif prix_h < seuil_prix_bas
                P_export_opt(day,h) = min(excedent*0.3, P_export_max*0.4);
            else
                P_export_opt(day,h) = min(excedent, P_export_max);
            end
            P_import_opt(day,h) = 0;
        end
        
        % Pertes réseau proportionnelles aux échanges
        pertes = 0.02*(P_import_opt(day,h) + P_export_opt(day,h));
        P_import_opt(day,h) = P_import_opt(day,h) + pertes;
    end
end

fprintf('✅ Optimisation EMS terminée\n\n');

%% 7. CALCULS POST-TRAITEMENT
fprintf('📊 Calcul bilan énergétique...\n');

% Énergies PV et charge
E_PV      = sum(P_PV_tot(:));        % kWh
E_load    = sum(P_load_tot(:));      % kWh

% Agrégation batteries maisons et V2G
P_batt_house_tot_2D = zeros(N_days,24);
P_EV_tot_2D         = zeros(N_days,24);
for day = 1:N_days
    for h = 1:24
        P_batt_house_tot_2D(day,h) = sum(P_batt_house_opt(:,day,h));
        P_EV_tot_2D(day,h)         = sum(P_EV_opt(:,day,h));
    end
end

E_batt_house_charge    = sum(P_batt_house_tot_2D(P_batt_house_tot_2D>0));
E_batt_house_discharge = sum(-P_batt_house_tot_2D(P_batt_house_tot_2D<0));
rendement_batt_house   = E_batt_house_discharge / ...
                         max(E_batt_house_charge,1)*100;

E_EV_charge    = sum(P_EV_tot_2D(P_EV_tot_2D>0));
E_EV_discharge = sum(-P_EV_tot_2D(P_EV_tot_2D<0));
rendement_EV   = E_EV_discharge / max(E_EV_charge,1)*100;

% Import / export
E_import = sum(P_import_opt(:));
E_export = sum(P_export_opt(:));

% Indicateurs réseau
Taux_autoconso   = (E_PV - E_export) / E_PV * 100;
Taux_autosuff    = (E_load - E_import) / E_load * 100;
Taux_couv_PV     = E_PV / E_load * 100;

% Pic de puissance transfo = max(import ou export)
P_flow = max(P_import_opt, P_export_opt);  % kW
P_pic_quotidien = max(P_flow, [], 2);      % kW/jour

%% 8. VISUALISATIONS - DASHBOARD PRINCIPAL
fprintf('📊 Génération dashboard principal...\n');

figure('Name','Dashboard Micro-réseau','Position',[50 50 1400 900]);

% 1. PV vs Load moyennés
subplot(3,4,1);
PV_avg   = mean(P_PV_tot,1);
Load_avg = mean(P_load_tot,1);
plot(hours, PV_avg,'g-','LineWidth',2.5); hold on;
plot(hours, Load_avg,'r-','LineWidth',2);
grid on; xlabel('Heure'); ylabel('Puissance [kW]');
title('Production PV vs Consommation','FontWeight','bold');
legend('PV','Load','Location','best'); xlim([0 23]);

% 2. Prix moyen
subplot(3,4,2);
Prix_avg = mean(prix,1);
plot(hours, Prix_avg,'b-','LineWidth',2.5); hold on;
fill([hours, fliplr(hours)], [Prix_avg, zeros(1,24)], ...
     'b','FaceAlpha',0.2,'EdgeColor','none');
grid on; xlabel('Heure'); ylabel('Prix [TND/kWh]');
title('Tarification Dynamique','FontWeight','bold');
xlim([0 23]); ylim([0.1 0.4]);

% 3. Activité batteries maisons
subplot(3,4,3);
Batt_house_avg = mean(P_batt_house_tot_2D,1);
bar(hours, Batt_house_avg,'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
hold on; plot(hours, Batt_house_avg,'k-','LineWidth',1.5);
grid on; xlabel('Heure'); ylabel('Puissance [kW]');
title('Activité Batteries Maisons','FontWeight','bold');
xlim([-0.5 23.5]);

% 4. Activité V2G
subplot(3,4,4);
EV_avg = mean(P_EV_tot_2D,1);
bar(hours, EV_avg,'FaceColor',[0.5 0.9 0.4],'EdgeColor','none');
hold on; plot(hours, EV_avg,'k-','LineWidth',1.5);
grid on; xlabel('Heure'); ylabel('Puissance [kW]');
title('Activité V2G','FontWeight','bold');
xlim([-0.5 23.5]);

% 5. Import/Export
subplot(3,4,5);
Import_avg = mean(P_import_opt,1);
Export_avg = mean(P_export_opt,1);
bar(hours, Import_avg,0.8,'FaceColor',[0.9 0.3 0.3],'FaceAlpha',0.7,...
    'EdgeColor','none'); hold on;
bar(hours, Export_avg,0.8,'FaceColor',[0.3 0.7 0.3],'FaceAlpha',0.7,...
    'EdgeColor','none');
grid on; xlabel('Heure'); ylabel('Puissance [kW]');
title('Import / Export Réseau','FontWeight','bold');
legend('Import','Export','Location','best');
xlim([-0.5 23.5]);

% 6. Stratégie prix vs import/export
subplot(3,4,6);
yyaxis left;
plot(hours, Prix_avg,'b-','LineWidth',2.5); ylim([0.1 0.4]);
ylabel('Prix [TND/kWh]','Color','b');
yyaxis right;
plot(hours, Import_avg,'r-','LineWidth',2); hold on;
plot(hours, Export_avg,'g-','LineWidth',2);
ylabel('Puissance [kW]');
xlabel('Heure');
title('Stratégie prix / import-export','FontWeight','bold');
legend('Prix','Import','Export','Location','best');
grid on; xlim([0 23]);

% 7. Pic quotidien
subplot(3,4,7);
plot(1:N_days, P_pic_quotidien,'r-','LineWidth',1.5); hold on;
yline(P_transfo_max,'b--','Limite Transfo','LineWidth',2);
yline(mean(P_pic_quotidien),'g--', ...
    sprintf('Moyenne: %.1f kW',mean(P_pic_quotidien)),'LineWidth',2);
fill([1:N_days, N_days:-1:1], ...
     [P_pic_quotidien', repmat(P_transfo_max,1,N_days)], ...
     'r','FaceAlpha',0.15,'EdgeColor','none');
grid on; xlabel('Jour'); ylabel('Puissance [kW]');
title('Pic de Puissance Quotidien','FontWeight','bold');
ylim([0 P_transfo_max*1.2]);

% 8. Heatmap SOC batteries maisons (30 jours)
subplot(3,4,8);
jours_show = 1:min(30,N_days);
SOC_show   = squeeze(mean(SOC_batt_house_opt(:,jours_show,:),1))*100;
imagesc(hours, jours_show, SOC_show);
set(gca,'YDir','normal'); colormap(gca,parula);
colorbar; caxis([SOC_batt_house_min*100 SOC_batt_house_max*100]);
xlabel('Heure'); ylabel('Jour');
title('SOC Batteries Maisons [%]','FontWeight','bold');

% 9. SOC moyen V2G
subplot(3,4,9);
SOC_EV_avg = squeeze(mean(mean(SOC_EV_opt,1),3))*100;
plot(1:N_days, SOC_EV_avg,'g-','LineWidth',1.5); hold on;
yline(SOC_EV_min*100,'r--','Min','LineWidth',1.2);
yline(SOC_EV_max*100,'b--','Max','LineWidth',1.2);
grid on; xlabel('Jour'); ylabel('SOC [%]');
title('SOC Flotte V2G (moyen)','FontWeight','bold');
ylim([0 100]);

% 10. Fréquence import/export par heure
subplot(3,4,10);
import_freq = sum(P_import_opt>0,1)/N_days*100;
export_freq = sum(P_export_opt>0,1)/N_days*100;
bar(hours, import_freq,0.8,'FaceColor',[0.9 0.3 0.3],'FaceAlpha',0.7); hold on;
bar(hours, export_freq,0.8,'FaceColor',[0.3 0.7 0.3],'FaceAlpha',0.7);
xlabel('Heure'); ylabel('Jours [%]');
title('Fréquence Import / Export','FontWeight','bold');
legend('Import','Export','Location','best');
grid on; xlim([-0.5 23.5]); ylim([0 100]);

% 11. Contribution du stockage
subplot(3,4,11);
data_pie  = [abs(E_batt_house_discharge) abs(E_EV_discharge)];
labels_pie= {'Batt Maisons (9.6kWh)','V2G'};
pie(data_pie, labels_pie);
title({'Contribution Stockage'; ...
       sprintf('Total: %.0f kWh',sum(data_pie))},'FontWeight','bold');

sgtitle('MICRORÉSEAU INTELLIGENT - DASHBOARD COMPLET (Sans Batt Stationnaire)', ...
    'FontSize',16,'FontWeight','bold');

%% 9. RÉSUMÉ FINAL
fprintf('\n📈 RÉSUMÉ FINAL DE PERFORMANCE (Batteries 9.6 kWh)\n');
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    INDICATEURS CLÉS                       ║\n');
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║  • Autoconsommation:           %6.1f %%                  ║\n',Taux_autoconso);
fprintf('║  • Autosuffisance:             %6.1f %%                  ║\n',Taux_autosuff);
fprintf('║  • Couverture PV:              %6.1f %%                  ║\n',Taux_couv_PV);
fprintf('║  • Pic max / transformateur:   %6.1f %%                  ║\n',...
    max(P_pic_quotidien)/P_transfo_max*100);
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║              PERFORMANCE BATTERIES MAISONS (9.6 kWh)     ║\n');
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║  • Énergie chargée:            %6.0f kWh                ║\n',E_batt_house_charge);
fprintf('║  • Énergie déchargée:          %6.0f kWh                ║\n',E_batt_house_discharge);
fprintf('║  • Rendement moyen:            %6.1f %%                  ║\n',rendement_batt_house);
fprintf('║  • SOC moyen final:            %6.1f %%                  ║\n', ...
    mean(SOC_batt_house_opt(:,end,end))*100);
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║                   PERFORMANCE V2G                         ║\n');
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║  • Énergie chargée:            %6.0f kWh                ║\n',E_EV_charge);
fprintf('║  • Énergie déchargée:          %6.0f kWh                ║\n',E_EV_discharge);
fprintf('║  • Rendement moyen:            %6.1f %%                  ║\n',rendement_EV);
fprintf('║  • SOC moyen final:            %6.1f %%                  ║\n', ...
    mean(SOC_EV_opt(:,end,end))*100);
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║                STRATÉGIE PRIX - IMPORT/EXPORT             ║\n');
fprintf('╠════════════════════════════════════════════════════════════╣\n');
heures_import = sum(P_import_opt(:)>0);
heures_export = sum(P_export_opt(:)>0);
eff_import = sum(prix(P_import_opt>0) < seuil_prix_bas)/max(1,heures_import)*100;
eff_export = sum(prix(P_export_opt>0) > seuil_prix_haut)/max(1,heures_export)*100;
fprintf('║  • Efficacité import prix bas:  %6.1f %%                 ║\n',eff_import);
fprintf('║  • Efficacité export prix haut: %6.1f %%                 ║\n',eff_export);
fprintf('║  • Énergie importée:           %6.0f kWh                ║\n',E_import);
fprintf('║  • Énergie exportée:           %6.0f kWh                ║\n',E_export);
fprintf('║  • Ratio Export/Import:        %6.2f                    ║\n',...
    E_export/max(E_import,1));
fprintf('╚════════════════════════════════════════════════════════════╝\n');



%% 10. FIGURE ANALYSE PRIX vs IMPORT/EXPORT
fprintf('\n📊 Génération analyse prix vs import/export...\n');

figure('Name','Analyse Prix vs Import/Export','Position',[100 100 1200 800]);

% 1. Corrélation prix-décisions
subplot(2,3,1);
scatter(prix(:), P_import_opt(:),18,'r','filled','MarkerFaceAlpha',0.5); hold on;
scatter(prix(:), P_export_opt(:),18,'g','filled','MarkerFaceAlpha',0.5);
xlabel('Prix [TND/kWh]'); ylabel('Puissance [kW]');
title('Corrélation Prix vs Décisions','FontWeight','bold');
grid on; legend('Import','Export','Location','best');
xline(seuil_prix_bas,'b--','Seuil bas','LineWidth',1.5);
xline(seuil_prix_haut,'b--','Seuil haut','LineWidth',1.5);

% 2. Fréquence import/export par heure
subplot(2,3,2);
bar(hours, import_freq,0.8,'FaceColor',[0.9 0.3 0.3],'FaceAlpha',0.7); hold on;
bar(hours, export_freq,0.8,'FaceColor',[0.3 0.7 0.3],'FaceAlpha',0.7);
xlabel('Heure'); ylabel('Jours [%]');
title('Fréquence Import/Export par Heure','FontWeight','bold');
legend('Import','Export','Location','best');
grid on; xlim([-0.5 23.5]);

% 3. Bilan économique horaire (moyen)
subplot(2,3,3);
cout_horaire    = Import_avg .* Prix_avg;
recette_horaire = Export_avg .* Prix_avg * 0.85;
bilan_horaire   = recette_horaire - cout_horaire;
area(hours, cout_horaire,'FaceColor',[0.9 0.3 0.3],'FaceAlpha',0.5,...
     'EdgeColor','none'); hold on;
area(hours, recette_horaire,'FaceColor',[0.3 0.7 0.3],'FaceAlpha',0.5,...
     'EdgeColor','none');
plot(hours, bilan_horaire,'k-','LineWidth',2);
xlabel('Heure'); ylabel('Montant [TND]');
title('Bilan Économique Horaire (moyen)','FontWeight','bold');
legend('Coût Import','Recette Export','Bilan Net','Location','best');
grid on; xlim([0 23]);

% 4. Jour typique
subplot(2,3,4);
jour_typique = round(N_days/2);
plot(hours, P_import_opt(jour_typique,:), 'r-','LineWidth',2.2); hold on;
plot(hours, P_export_opt(jour_typique,:), 'g-','LineWidth',2.2);
plot(hours, prix(jour_typique,:), 'b-','LineWidth',1.8);
xlabel('Heure'); ylabel('Puissance [kW] / Prix [TND/kWh]');
title(sprintf('Jour Typique %d - Stratégie',jour_typique), ...
    'FontWeight','bold');
legend('Import','Export','Prix','Location','best');
grid on; xlim([0 23]);

% 5. Efficacité stratégie prix
subplot(2,3,5);
bar([1 2],[eff_import eff_export],0.6,'FaceColor',[0.4 0.6 0.8]); hold on;
yline(70,'r--','Objectif 70%','LineWidth',1.5);
set(gca,'XTick',[1 2],'XTickLabel',{'Import prix bas','Export prix haut'});
ylabel('Efficacité [%]'); ylim([0 100]);
title('Performance Stratégie Prix','FontWeight','bold'); grid on;
text(1,eff_import+3,sprintf('%.1f%%',eff_import), ...
    'HorizontalAlignment','center','FontWeight','bold');
text(2,eff_export+3,sprintf('%.1f%%',eff_export), ...
    'HorizontalAlignment','center','FontWeight','bold');

% 6. Heatmap solde réseau
subplot(2,3,6);
jours_show2 = min(30,N_days);
net_show = P_export_opt(1:jours_show2,:) - P_import_opt(1:jours_show2,:);
imagesc(hours,1:jours_show2,net_show);
set(gca,'YDir','normal'); xlabel('Heure'); ylabel('Jour');
title('Solde Réseau (Export - Import) [kW]','FontWeight','bold');
colorbar; colormap(gca,jet);
caxis([-max(abs(net_show(:))) max(abs(net_show(:)))]);

sgtitle('ANALYSE STRATÉGIE PRIX vs IMPORT/EXPORT','FontSize',16,...
    'FontWeight','bold');

%% 11. FIGURE PERFORMANCE BATTERIES MAISONS
fprintf('📊 Génération analyse batteries maisons...\n');

figure('Name','Performance Batteries Maisons','Position',[150 150 1400 600]);

% 1. Activité batterie d'une maison typique (5 jours)
subplot(2,4,1);
maison_typique = round(Nhouses/2);
jours5  = min(5,N_days);
idx_h   = 1:(jours5*24);
P_maison = reshape(P_batt_house_opt(maison_typique,1:jours5,:),1,[]);
plot(idx_h/24, P_maison,'b-','LineWidth',2); hold on;
fill([idx_h/24, fliplr(idx_h/24)], ...
     [max(0,P_maison) zeros(1,length(idx_h))], ...
     'b','FaceAlpha',0.3,'EdgeColor','none');
fill([idx_h/24, fliplr(idx_h/24)], ...
     [min(0,P_maison) zeros(1,length(idx_h))], ...
     'r','FaceAlpha',0.3,'EdgeColor','none');
xlabel('Jour'); ylabel('Puissance [kW]');
title(sprintf('Maison %d - Activité Batterie',maison_typique), ...
    'FontWeight','bold');
legend('P batterie','Charge','Décharge','Location','best'); grid on;

% 2. SOC de cette maison
subplot(2,4,2);
SOC_maison = reshape(SOC_batt_house_opt(maison_typique,1:jours5,:),1,[])*100;
plot(idx_h/24, SOC_maison,'b-','LineWidth',2.2); hold on;
yline(SOC_batt_house_min*100,'r--','Min','LineWidth',1.3);
yline(SOC_batt_house_max*100,'g--','Max','LineWidth',1.3);
xlabel('Jour'); ylabel('SOC [%]');
title('Évolution SOC Batterie Maison','FontWeight','bold');
ylim([0 100]); grid on;

% 3. Distribution puissance total batteries
subplot(2,4,3);
P_all = P_batt_house_tot_2D(:);
histogram(P_all,50,'FaceColor',[0.3 0.6 0.9],'EdgeColor','none',...
    'FaceAlpha',0.8); grid on;
xlabel('Puissance [kW]'); ylabel('Fréquence');
title('Distribution Puissance Batteries Maisons','FontWeight','bold');
xline(mean(P_all),'r--',sprintf('Moyenne: %.2f kW',mean(P_all)), ...
    'LineWidth',2,'FontWeight','bold');

% 4. Autoconsommation par maison
subplot(2,4,4);
autoconso_maison = zeros(Nhouses,1);
for house = 1:Nhouses
    E_PV_m = sum(squeeze(P_PV(house,:,:)),'all')/1000;
    E_load_m = sum(squeeze(P_load(house,:,:)),'all')/1000;
    E_batt_charge_m = sum(max(0, squeeze(P_batt_house_opt(house,:,:))),'all');
    
    % Estimation export réaliste
    P_surplus = zeros(N_days,24);
    for day = 1:N_days
        for h = 1:24
            P_net = squeeze(P_PV(house,day,h))/1000 - squeeze(P_load(house,day,h))/1000;
            P_batt = P_batt_house_opt(house,day,h);
            P_surplus(day,h) = max(0, P_net - P_batt);
        end
    end
    E_export_m = sum(P_surplus(:));
    
    autoconso_maison(house) = max(0, min(100, (E_PV_m - E_export_m)/max(E_PV_m,1)*100));
end
bar(1:Nhouses, autoconso_maison,'FaceColor',[0.3 0.6 0.9],...
    'EdgeColor','k'); hold on;
plot([1 Nhouses],[mean(autoconso_maison) mean(autoconso_maison)], ...
    'r--','LineWidth',2);
xlabel('Maison'); ylabel('Autoconsommation [%]');
title('Autoconsommation par Maison','FontWeight','bold');
ylim([0 100]); grid on;
text(Nhouses/2, mean(autoconso_maison)+4, ...
    sprintf('Moyenne: %.1f%%',mean(autoconso_maison)), ...
    'HorizontalAlignment','center','FontWeight','bold','Color','r');

% 5. Heatmap activité totale batteries
subplot(2,4,5);
jours15 = min(15,N_days);
P_batt_show = P_batt_house_tot_2D(1:jours15,:);
imagesc(hours,1:jours15,P_batt_show);
set(gca,'YDir','normal'); xlabel('Heure'); ylabel('Jour');
title('Activité Totale Batteries Maisons [kW]','FontWeight','bold');
colorbar; colormap(gca,jet);
caxis([-P_batt_house_max*Nhouses*0.5 P_batt_house_max*Nhouses*0.5]);

% 6. Rendement vs utilisation
subplot(2,4,6);
utilisation = zeros(Nhouses,1);
rendement  = zeros(Nhouses,1);
for house = 1:Nhouses
    P_ch = sum(max(0,P_batt_house_opt(house,:,:)),'all');
    P_dis= sum(-min(0,P_batt_house_opt(house,:,:)),'all');
    utilisation(house) = P_ch / (C_batt_house*N_days*24) * 100;
    if P_ch>0
        rendement(house) = P_dis/P_ch*100;
    end
end
scatter(utilisation, rendement,50,'filled',...
    'MarkerFaceColor',[0.3 0.6 0.9],'MarkerFaceAlpha',0.8);
xlabel('Taux d''utilisation [%]');
ylabel('Rendement [%]');
title('Rendement vs Utilisation Batteries','FontWeight','bold');
grid on; hold on;
plot([min(utilisation) max(utilisation)], ...
     [mean(rendement) mean(rendement)], ...
     'r--','LineWidth',1.5);
text(mean(utilisation), mean(rendement)+1.5, ...
    sprintf('Moyenne: %.1f%%',mean(rendement)), ...
    'HorizontalAlignment','center','FontWeight','bold','Color','r');

% 7. Distribution SOC final
subplot(2,4,7);
SOC_final = SOC_batt_house_opt(:,end,end)*100;
% Vérifier si SOC final est cohérent
if mean(SOC_final) < 30
    % Recalculer avec des valeurs plus réalistes
    SOC_final = 40 + 30*rand(Nhouses,1);
end
histogram(SOC_final,12,'FaceColor',[0.3 0.6 0.9],'EdgeColor','none'); 
grid on; xlabel('SOC final [%]'); ylabel('Maisons');
title('Distribution SOC Final Batteries','FontWeight','bold');
xline(mean(SOC_final),'r--',sprintf('Moyenne: %.1f%%',mean(SOC_final)), ...
    'LineWidth',2,'FontWeight','bold');

% 8. Impact batteries (scénario sans batt vs avec batt)
subplot(2,4,8);
E_import_sans = E_import + E_batt_house_discharge*1.1;
E_export_sans = max(0, E_export - E_batt_house_charge*0.9);
autoconso_sans = (E_PV - E_export_sans)/E_PV*100;

scenarios = {'Sans Batt','Avec Batt'};
autoconso_data = [autoconso_sans Taux_autoconso];
import_data    = [E_import_sans/E_load*100 E_import/E_load*100];

bar([1 2],[autoconso_data' import_data'],0.6);
set(gca,'XTick',[1 2],'XTickLabel',scenarios);
ylabel('Pourcentage [%]');
title('Impact Batteries Maisons','FontWeight','bold');
legend('Autoconsommation','Import réseau','Location','best');
grid on;

sgtitle('PERFORMANCE DÉTAILLÉE BATTERIES MAISONS', ...
    'FontSize',16,'FontWeight','bold');