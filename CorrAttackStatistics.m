clear all; close all; clc; % limpiamos el Entorno de trabajo
%%

Vaults= dir(strcat(pwd,'\ExpOctubre\Vaults\*.txt'));

cont = 0;
huf = 1;
for hu=1:length(Vaults)
cont = cont + 1;   
if cont == 9
    huf = huf+1;
    cont = 1;
end
im_no(huf,cont) = string(Vaults(hu).name(6:10));
inds(huf,cont) = cont;
end

[s1_ino,s2_ino] = size(im_no);
for i =1:s1_ino
    combs(i,:, :) =  nchoosek(inds(i,:),2);
end
[s1,s2,s3] = size(combs);

for r = 1:s1
   for t = 1:s2
       comb_names(r,t,:) = im_no(r,combs(r,t,:));
   end
end


%%
distancias=linspace(2000,6000,4);
% distancias = [1800];

for d = 1:length(distancias)
    fprintf('\n\n\t\t\tVoy en el experimento: %d \n\n', d);
    [scn1,scn2,scn3]=size(comb_names);
    m_n = 0;
    cont_combs = 0;
    conts = 0;
    for vu = 1:scn1
        for vc = 1:scn2
            conts = conts + 1;
            
            VaultFijName = comb_names(vu,vc,1);
            VaultMovName = comb_names(vu,vc,2);

            VaultFij = load(strcat('ExpOctubre\Vaults\Vault',VaultFijName, '.txt'));  %%%% Vault Extraída de Sistema 1
            VaultMov = load(strcat('ExpOctubre\Vaults\Vault',VaultMovName, '.txt'));  %%%% Vault Extraída de Sistema 2



            Real_XY1 = load(strcat('ExpOctubre\Real_XYs\Real_XY', VaultFijName, '.txt'));  % Puntos genuinos VaultFij
            Real_XY2 = load(strcat('ExpOctubre\Real_XYs\Real_XY', VaultMovName, '.txt'));  % Puntos genuinos VaultMov
            Chaff_Data1 = load(strcat('ExpOctubre\Chaff_Datas\Chaff_Data', VaultFijName, '.txt'));  % Puntos de ruido VaultFij
            Chaff_Data2 = load(strcat('ExpOctubre\Chaff_Datas\Chaff_Data', VaultMovName,'.txt'));  % Puntos de ruido VaultMov



            min_dist = distancias(d);%2500;%1800;%2000;%2250; % distancia minima para identificar puntos coincidentes
            tran_dist = 500;%Distancias a recorrer para buscar un mejor alineamiento
            min_combs = 100000;
            
            VaultMov = VaultMov'; 
            VaultXt = VaultMov(1,:) - mean(VaultMov(1,:));
            VaultYt = VaultMov(2,:) - mean(VaultMov(2,:));
            VaultT = [VaultXt;VaultYt];

            rotattempts = 1500; % Se define el número de rotaciones. 
            rots = linspace(0,2,rotattempts); % creamos intervalos para rotar  tomando 


            transattempts = 1000; % número de translaciones que se harán. 

            for i = 1:rotattempts % Se ejecutan las rotaciones

                R = [cos(3.1416*rots(i)) -sin(3.1416*rots(i)); sin(3.1416*rots(i)) cos(3.1416*rots(i))]; 
                DataRot = R * VaultT; 
                VaultMovRot = [DataRot(1,:)+ mean(VaultMov(1,:));DataRot(2,:)+ mean(VaultMov(2,:))]; 
                VaultMovRot = VaultMovRot'; 

                [Vault_indiceV1,Real_PointsV1,Vault_indiceV2,Real_PointsV2]=distancia_VaultsCA(VaultFij,VaultMovRot, min_dist,0);

                MatchPoints(i) = length(Vault_indiceV1);
            end


            maxmp = max(MatchPoints); % se saca el número mayor de puntos coincidentes entre todas las rotaciones
            maxmpidx = find(MatchPoints==maxmp); % Se saca que rotaciones obtuvieron el maxmp 

            for j = 1:length(maxmpidx)
            %     fprintf('\n Iteración j: %d  ', j)
                % Se rota la bóveda al j un grado  dentro de  maxmpidx
                R = [cos(3.1416*rots(maxmpidx(j))) -sin(3.1416*rots(maxmpidx(j))); sin(3.1416*rots(maxmpidx(j))) cos(3.1416*rots(maxmpidx(j)))];
                DataRot = R * VaultT;
                VaultMovRot = [DataRot(1,:)+ mean(VaultMov(1,:));DataRot(2,:)+ mean(VaultMov(2,:))];
                VaultMovRot = VaultMovRot';
                
                % Se harán para cada rotación transattempts intentos de translación
                % para buscar el mayor numero de puntos coincidentes
                for k = 1:transattempts
            %         fprintf('\n Iteración k: %d  ', k)
                    % se calculan los puntos coincidentes iniciales cada que entra al
                    % for ya que si mejoran en cierta interación se guardan en la
                    % vóbeda que entra aqui, ya con la translación de mejora hecha. Si
                    % no mejora, se mantiene igual. 
                    [Vault_indiceV1,Real_PointsV1,Vault_indiceV2,Real_PointsV2]=distancia_VaultsCA(VaultFij,VaultMovRot, min_dist,0);
                    
                    % Se trnaladan los funtos de VaultMov ya rotada en las 8
                    % direcciones posibles para los pixeles. 
                    dir1 = [VaultMovRot(:,1)+tran_dist,VaultMovRot(:,2)]; % Derecha
                    dir2 = [VaultMovRot(:,1)-tran_dist,VaultMovRot(:,2)]; % Izquierda
                    dir3 = [VaultMovRot(:,1),VaultMovRot(:,2)+tran_dist]; % Arriba
                    dir4 = [VaultMovRot(:,1),VaultMovRot(:,2)-tran_dist]; % Abajo
                    dir5 = [VaultMovRot(:,1)+tran_dist,VaultMovRot(:,2)+tran_dist]; % Arriba - Derecha
                    dir6 = [VaultMovRot(:,1)+tran_dist,VaultMovRot(:,2)-tran_dist]; % Abajo - Derecha
                    dir7 = [VaultMovRot(:,1)-tran_dist,VaultMovRot(:,2)+tran_dist]; % Arriba - Izquierda
                    dir8 = [VaultMovRot(:,1)-tran_dist,VaultMovRot(:,2)-tran_dist]; % Abajo - Izquierda
                    dirs = cat(3,dir1,dir2,dir3,dir4,dir5,dir6,dir7,dir8); % Se guardan las translaciones hechas
                    [s1,s2,s3] = size(dirs);
                    
                    % Ahora se calculará la distacia para cada translación y se
                    % obtendra la mejor, si es que hay. La que de mejores numeros de
                    % coincidencias
                    mk = 0; 
                    
                    for h = 1:length(s3)
                        [Vault_indiceV1T,Real_PointsV1T,Vault_indiceV2T,Real_PointsV2T]=distancia_VaultsCA(VaultFij,dirs(:,:,h), min_dist,0);
                        % Se verifica si para la h traslación para la k rotación da mas
                        % puntos coincidentes y si si, se establece como la marca a
                        % superar y se guarda la bóveda en esas condiciones. 
                        if length(Vault_indiceV1T) > mk 
                           mk =  length(Vault_indiceV1T);
                           best_trans = dirs(:,:,h);
                        else
                            best_trans = [];
                        end
                    end
                    
                    if isempty(best_trans)
                        continue
                    else  % se checa si si mejoro la coincidencia con alguna rotación
                        [Vault_indiceV1BT,Real_PointsV1BT,Vault_indiceV2BT,Real_PointsV2BT]=distancia_VaultsCA(VaultFij,best_trans, min_dist,0);
                        % se verifica que realmente sea mejor a la del inicio del for para la k rotación. 
                        if length(Vault_indiceV1BT) > length(Vault_indiceV1)
                            % Si si, se guarda como la bóveda original para continuar con el for y ver si se puede mejorar aun mas
                            VaultMovRot = best_trans; 
                        end
                    end
                end
                clear b_trans4rot
                b_trans4rot(:,:, j) = VaultMovRot; % Se guarda cada rotación de las mejores con translaciones hasta que mas coincidencias diera
            end
            
            % Ahora de las rotaciones con tranlaciones hasta el mejor punto, se toma
            % solo la que mas puntos coincidentes tenga de entre estas.
            [sb1,sb2,sb3] = size(b_trans4rot);
            mk2 = 0;
            for f =1:sb3
                [Vault_indiceV1BT4R,Real_PointsV1BT4R,Vault_indiceV2BT4R,Real_PointsV2BT4R]=distancia_VaultsCA(VaultFij,b_trans4rot(:,:,f), min_dist,0);
                if length(Vault_indiceV1BT4R)> mk2
                   mk2 = length(Vault_indiceV1BT4R);
                   VaultMovRotf = b_trans4rot(:,:,f); % Se establece entonces la mejor VaultMov con transformaciones obtenida
                end
            end

%             mk2 = 0;
%             for f =1:length(maxmpidx)
%                 R = [cos(3.1416*rots(maxmpidx(f))) -sin(3.1416*rots(maxmpidx(f))); sin(3.1416*rots(maxmpidx(f))) cos(3.1416*rots(maxmpidx(f)))];
%                 DataRot = R * VaultT;
%                 VaultMovRot = [DataRot(1,:)+ mean(VaultMov(1,:));DataRot(2,:)+ mean(VaultMov(2,:))];
%                 VaultMovRot = VaultMovRot';
% 
%                 [Vault_indiceV1BT4R,Real_PointsV1BT4R,Vault_indiceV2BT4R,Real_PointsV2BT4R]=distancia_VaultsCA(VaultFij,VaultMovRot, min_dist,0);
%                 if length(Vault_indiceV1BT4R)> mk2
%                    mk2 = length(Vault_indiceV1BT4R);
%                    VaultMovRotf = VaultMovRot; % Se establece entonces la mejor VaultMov con transformaciones obtenida
%                 end
%             end
%             VaultMovRot = VaultMovRotf;

            [Vault_indiceV1Check,Real_PointsV1Check,Vault_indiceV2Check,Real_PointsV2Check]=distancia_VaultsCA(VaultFij,VaultMovRotf, min_dist,0);
            g_enc = 0;

            for c = 1:length(Real_PointsV1Check)
                for c2 = 1:length(Real_XY1)
                    if Real_PointsV1Check(c,1) == Real_XY1(c2,1) && Real_PointsV1Check(c,2) == Real_XY1(c2,2)
                        g_enc = g_enc + 1;
                    end
                end
            end

            if g_enc >= 9
                m_n = m_n + 1;
                ca_ex(m_n,1) = VaultFijName;
                ca_ex(m_n,2) = VaultMovName;
                combs_pos = factorial(length(Real_PointsV1Check))/(factorial(9)*factorial(length(Real_PointsV1Check)-9));
                combs_prom = combs_pos / (factorial(g_enc)/(factorial(9)*factorial(g_enc-9)));
                fprintf('Para %s y %s las combinaciones posibles son: %d y las promedio son %d \n', VaultFijName,VaultMovName,combs_pos,combs_prom);
                if combs_prom < min_combs
                    cont_combs =  cont_combs +1;
                    combs_ft(cont_combs,:) = [VaultFijName, VaultMovName, combs_pos, combs_prom];
%                     fprintf('Para %s y %s las combinacione sposibles son: %d y las promedio son %d \n', VaultFijName,VaultMovName,combs_pos,combs_prom);
                    fprintf('#------------- Combinación de huellas viable encontrada ------------------# \n')
                end
            end
            fprintf('Los puntos genuinos entre %s y %s encontrados son: %d de entre %d en total coincidentes. \n', VaultFijName,VaultMovName,g_enc,length(Real_PointsV1Check));
            
            resultados(d,conts,:) = [g_enc;length(Real_PointsV1Check)];
        end 
    end

%        fprintf('Voy en el experimento: %d \n', d);
      Tm_n(d) = m_n;
%     fprintf('El total de casos que tendría exito es: %d \n', m_n);

%     fprintf('las bóvedas que al compararse tendria exito el ataque son: \n');

%     display(ca_ex)
end

save('ExpOctubre\Statistics\ResCorr4mindist1000trans.mat','resultados')
writematrix(resultados,'ExpOctubre\Statistics\ResCorr4mindist1000trans.csv')
% save('ExpOctubre\Statistics\CombsPos_mindist1800_SoloGiro.mat','combs_ft')

%%
[srp1,srp2,srp3] = size(resultados);
% for rp = 1:srp1
%     crp = 0;
%     for rp2 = 1:srp2
%        if resultados(rp,rp2,1) >= 9
%            crp = crp +1;
%            res_pos(rp,crp,1) = resultados(rp,rp2,1);
%            res_pos(rp,crp,2) = resultados(rp,rp2,2);
%        end
%     end
% end

crp = 0;
crp2 = 0;
crp3 = 0;
crp4 = 0;
crp5 = 0;
crp6 = 0;
crp7 = 0;
crp8 = 0;
crp9 = 0;
crp10 = 0;

for rp2 = 1:srp2
   if resultados(1,rp2,1) >= 9
       crp = crp +1;
       res_pos1(crp,1) = resultados(1,rp2,1);
       res_pos1(crp,2) = resultados(1,rp2,2);
   end
   if resultados(2,rp2,1) >= 9
       crp2 = crp2 +1;
       res_pos2(crp2,1) = resultados(2,rp2,1);
       res_pos2(crp2,2) = resultados(2,rp2,2);
   end
   if resultados(3,rp2,1) >= 9
       crp3 = crp3 +1;
       res_pos3(crp3,1) = resultados(3,rp2,1);
       res_pos3(crp3,2) = resultados(3,rp2,2);
   end
   if resultados(4,rp2,1) >= 9
       crp4 = crp4 +1;
       res_pos4(crp4,1) = resultados(4,rp2,1);
       res_pos4(crp4,2) = resultados(4,rp2,2);
   end
%    if resultados(5,rp2,1) >= 9
%        crp5 = crp5 +1;
%        res_pos5(crp5,1) = resultados(5,rp2,1);
%        res_pos5(crp5,2) = resultados(5,rp2,2);
%    end
%    if resultados(6,rp2,1) >= 9
%        crp6 = crp6 +1;
%        res_pos6(crp6,1) = resultados(6,rp2,1);
%        res_pos6(crp6,2) = resultados(6,rp2,2);
%    end
%    if resultados(7,rp2,1) >= 9
%        crp7 = crp7 +1;
%        res_pos7(crp7,1) = resultados(7,rp2,1);
%        res_pos7(crp7,2) = resultados(7,rp2,2);
%    end
%    if resultados(8,rp2,1) >= 9
%        crp8 = crp8 +1;
%        res_pos8(crp8,1) = resultados(8,rp2,1);
%        res_pos8(crp8,2) = resultados(8,rp2,2);
%    end
%    if resultados(9,rp2,1) >= 9
%        crp9 = crp9 +1;
%        res_pos9(crp9,1) = resultados(9,rp2,1);
%        res_pos9(crp9,2) = resultados(9,rp2,2);
%    end
%    if resultados(10,rp2,1) >= 9
%        crp10 = crp10 +1;
%        res_pos10(crp10,1) = resultados(10,rp2,1);
%        res_pos10(crp10,2) = resultados(10,rp2,2);
%    end
end

%%


% Av_G1 = 0;
% Av_N1 = 0;
% Av_A1=0;
Av_G1 = round(sum(res_pos1(:,1))/length(res_pos1(:,1)));
Av_N1 = round(sum(res_pos1(:,2))/length(res_pos1(:,2)));
Av_A1=(factorial(Av_N1)/(factorial(9)*factorial(Av_N1-9)))/(factorial(Av_G1)/(factorial(9)*factorial(Av_G1-9)));


Av_G2 = round(sum(res_pos2(:,1))/length(res_pos2(:,1)));
Av_N2 = round(sum(res_pos2(:,2))/length(res_pos2(:,2)));
Av_A2=(factorial(Av_N2)/(factorial(9)*factorial(Av_N2-9)))/(factorial(Av_G2)/(factorial(9)*factorial(Av_G2-9)));

Av_G3 = round(sum(res_pos3(:,1))/length(res_pos3(:,1)));
Av_N3 = round(sum(res_pos3(:,2))/length(res_pos3(:,2)));
Av_A3=(factorial(Av_N3)/(factorial(9)*factorial(Av_N3-9)))/(factorial(Av_G3)/(factorial(9)*factorial(Av_G3-9)));

Av_G4 = round(sum(res_pos4(:,1))/length(res_pos4(:,1)));
Av_N4 = round(sum(res_pos4(:,2))/length(res_pos4(:,2)));
Av_A4 =(factorial(Av_N4)/(factorial(9)*factorial(Av_N4-9)))/(factorial(Av_G4)/(factorial(9)*factorial(Av_G4-9)));

% Av_G5 = round(sum(res_pos5(:,1))/length(res_pos5(:,1)));
% Av_N5 = round(sum(res_pos5(:,2))/length(res_pos5(:,2)));
% Av_A5 =(factorial(Av_N5)/(factorial(9)*factorial(Av_N5-9)))/(factorial(Av_G5)/(factorial(9)*factorial(Av_G5-9)));
% 
% Av_G6 = round(sum(res_pos6(:,1))/length(res_pos6(:,1)));
% Av_N6 = round(sum(res_pos6(:,2))/length(res_pos6(:,2)));
% Av_A6 =(factorial(Av_N6)/(factorial(9)*factorial(Av_N6-9)))/(factorial(Av_G6)/(factorial(9)*factorial(Av_G6-9)));
% 
% Av_G7 = round(sum(res_pos7(:,1))/length(res_pos7(:,1)));
% Av_N7 = round(sum(res_pos7(:,2))/length(res_pos7(:,2)));
% Av_A7 = (factorial(Av_N7)/(factorial(9)*factorial(Av_N7-9)))/(factorial(Av_G7)/(factorial(9)*factorial(Av_G7-9)));
% 
% 
% Av_G8 = round(sum(res_pos8(:,1))/length(res_pos8(:,1)));
% Av_N8 = round(sum(res_pos8(:,2))/length(res_pos8(:,2)));
% Av_A8 =(factorial(Av_N8)/(factorial(9)*factorial(Av_N8-9)))/(factorial(Av_G8)/(factorial(9)*factorial(Av_G8-9)));
% 
% Av_G9 = round(sum(res_pos9(:,1))/length(res_pos9(:,1)));
% Av_N9 = round(sum(res_pos9(:,2))/length(res_pos9(:,2)));
% Av_A9 = (factorial(Av_N9)/(factorial(9)*factorial(Av_N9-9)))/(factorial(Av_G9)/(factorial(9)*factorial(Av_G9-9)));
% 
% Av_G10 = round(sum(res_pos10(:,1))/length(res_pos10(:,1)));
% Av_N10 = round(sum(res_pos10(:,2))/length(res_pos10(:,2)));
% Av_A10 = (factorial(Av_N10)/(factorial(9)*factorial(Av_N10-9)))/(factorial(Av_G10)/(factorial(9)*factorial(Av_G10-9)));

% for r = 1:10
%     
%     Av_G = round(sum(resultados(r,:,1))/length(resultados(r,:,1)));
%     Av_N = round(sum(resultados(r,:,2))/length(resultados(r,:,2)));
%     
%     if Av_G < 9
%         Av_A(r)= 0;
%     else
%         Av_A(r)=(factorial(Av_N)/(factorial(9)*factorial(Av_N-9)))/(factorial(Av_G)/(factorial(9)*factorial(Av_G-9)));
%     end
% end
%%

% TabRes = [distancias;[Av_N1,Av_N2,Av_N3,Av_N4,Av_N5,Av_N6,Av_N7,Av_N8,Av_N9,Av_N10];[Av_G1,Av_G2,Av_G3,Av_G4,Av_G5,Av_G6,Av_G7,Av_G8,Av_G9,Av_G10];Tm_n; [Av_A1,Av_A2,Av_A3,Av_A4,Av_A5,Av_A6,Av_A7,Av_A8,Av_A9,Av_A10]]
% TabRes = [distancias;[Av_N1,Av_N2,Av_N3,Av_N4,Av_N5];[Av_G1,Av_G2,Av_G3,Av_G4,Av_G5];Tm_n; [Av_A1,Av_A2,Av_A3,Av_A4,Av_A5]]
% TabRes = [distancias;[Av_N1,Av_N2,Av_N3];[Av_G1,Av_G2,Av_G3];Tm_n; [Av_A1,Av_A2,Av_A3]]
TabRes = [distancias;[Av_N1,Av_N2,Av_N3,Av_N4];[Av_G1,Av_G2,Av_G3,Av_G4];Tm_n; [Av_A1,Av_A2,Av_A3,Av_A4]]
save('ExpOctubre\Statistics\TabResCorr4mindist1000trans.mat','TabRes')
writematrix(TabRes,'ExpOctubre\Statistics\TabResCorr4mindist1000trans.csv')
%%

hFig1 = figure(1);
set(gcf,'position',get(0,'ScreenSize'))
hr = plot(VaultFij(:,1),VaultFij(:,2),'bv'), title('Puntos Reales Ambas Bovedas');
hold on;
plot(VaultMovRot(:,1), VaultMovRot(:,2),'rv');
%Aquí se ponen triangulos verdes sobre los azules genuinos que lograron ser identificados
%Por lo tanto los puntos de triangulo verde con una cruz 'X' negra encima
%son los puntos reales (genuinos) de VaultFij, que se lograron obtener con
%las rotaciones y tranlaciones
[Vault_indiceV1,Real_PointsV1,Vault_indiceV2,Real_PointsV2]=distancia_VaultsCA(VaultFij,VaultMovRot, min_dist,1);


plot(Real_XY1(:,1),Real_XY1(:,2),'kx'), title('9 Puntos Genuinos Identificados');
hold on;
plot(Real_XY2(:,1),Real_XY2(:,2),'rx');
plot(Chaff_Data1(:,1),Chaff_Data1(:,2),'ko');
plot(Chaff_Data2(:,1),Chaff_Data2(:,2),'ro');



    
  
