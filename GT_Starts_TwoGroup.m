%  Author: Dr. Rajanikant Panda (Email:bk.bme.rajanikant@gmail.com)
%  Description: This code used to compute the Normalised Graph
%  Mesures first for both the Control and Patient cohorts, followed by statistical analysis and visualization.

%% clear everything
clc; clear; close all

%% Group/Condition one

path = 'J:\fMRI\TINNITUS\Control\'
cd(path);
SUBJlist_Group1 = dir('dataanalysis*');

%% Absolute and Random and Normalised CC, PL, SW for Group1_ET data extraction
for i = 1:length(SUBJlist_Group1)
    %%
    SUBJname = SUBJlist_Group1(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group1_CC(i,:,:)= data.GT_clust_coeff;         
    Group1_PL(i,:,:)= data.GT_path_length;
    Group1_LE(i,:,:)= data.GT_local_eff;         
    Group1_GE(i,:,:)= data.GT_global_eff;
    Group1_Degree(i,:,:)= data.GT_degree;
    Group1_PC(i,:,:)= data.GT_participation_coeff;  
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group1_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group1_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group1_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group1_GE_rand(i,:,:,:)= data_rand.GT_global_eff ;  %GT_global_eff_rand;
    Group1_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group1_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
end

%%
Group1_CC_50=Group1_CC(:,1:20,:);                          
Group1_CC_rand_squ = squeeze(mean(Group1_CC_rand,3));        
AvgGroup1_CC=mean(mean(Group1_CC_50,3));

sparsity_CC_Group1_50 = (mean(Group1_CC_50,3));
sparsity_CC_rand_Group1 = mean(Group1_CC_rand_squ,3);
sparsity_CC_rand_Group1_50 = sparsity_CC_rand_Group1 (:,1:20);

Group1_PL_50=Group1_PL(:,:,1:20);
AvgGroup1_PL=squeeze(mean(Group1_PL_50));

Sparsity_PL_Group1_50=squeeze(Group1_PL_50);
Sparsity_PL_Group1_rand_squ = squeeze(mean(Group1_PL_rand,3)); 
Sparsity_PL_Group1_rand_50 = Sparsity_PL_Group1_rand_squ (:,1:20);
Group1_Degree_50=Group1_Degree(:,1:20,:);                             
AvgGroup1_Degree=mean(mean(Group1_Degree_50,3));
sparsity_Degree_Group1_50 = (mean(Group1_Degree_50,3));

Group1_PC_50=Group1_PC(:,1:20,:);                          
Group1_PC_rand_squ = squeeze(mean(Group1_PC_rand,3));        
AvgGroup1_PC=mean(mean(Group1_PC_50,3));

sparsity_PC_Group1_50 = (mean(Group1_PC_50,3));
sparsity_PC_rand_Group1 = mean(Group1_PC_rand_squ,3);
sparsity_PC_rand_Group1_50 = sparsity_PC_rand_Group1 (:,1:20);

%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            sparsity_CC_normalised_Group1(i,j) = sparsity_CC_Group1_50(i,j)/sparsity_CC_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            sparsity_PL_normalised_Group1(i,j) = Sparsity_PL_Group1_50(i,j)/Sparsity_PL_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            SmallWorldNess_Group1(i,j) = sparsity_CC_normalised_Group1(i,j)/sparsity_PL_normalised_Group1(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            sparsity_PC_normalised_Group1(i,j) = sparsity_PC_Group1_50(i,j)/sparsity_PC_rand_Group1_50(i,j); 
    end
end

%%
Group1_LE_50=Group1_LE(:,1:20,:);                          
Group1_LE_rand_squ = squeeze(mean(Group1_LE_rand,3));        
AvgGroup1_LE=mean(mean(Group1_LE_50,3));

sparsity_LE_Group1_50 = (mean(Group1_LE_50,3));
sparsity_LE_rand_Group1 = mean(Group1_LE_rand_squ,3);
sparsity_LE_rand_Group1_50 = sparsity_LE_rand_Group1 (:,1:20);

Group1_GE_50=Group1_GE(:,1:20);
AvgGroup1_GE=squeeze(mean(Group1_GE_50));

Sparsity_GE_Group1_50=squeeze(Group1_GE_50);
Sparsity_GE_Group1_rand_squ = squeeze(mean(Group1_GE_rand,3)); 
Sparsity_GE_Group1_rand_50 = Sparsity_GE_Group1_rand_squ (:,1:20); 

%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            sparsity_LE_normalised_Group1(i,j) = sparsity_LE_Group1_50(i,j)/sparsity_LE_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:20; 
            sparsity_GE_normalised_Group1(i,j) = Sparsity_GE_Group1_50(i,j)/Sparsity_GE_Group1_rand_50(i,j); 
    end
end
%% Group/Condition two study

path = 'J:\fMRI\TINNITUS\Patients\'
cd(path);
SUBJlist_Group2 = dir('dataanalysis*');
%%
for i = 1:length(SUBJlist_Group2)
    SUBJname = SUBJlist_Group2(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group2_CC(i,:,:)= data.GT_clust_coeff;         
    Group2_PL(i,:,:)= data.GT_path_length;
    Group2_LE(i,:,:)= data.GT_local_eff;         
    Group2_GE(i,:,:)= data.GT_global_eff;
    Group2_Degree(i,:,:)= data.GT_degree;
    Group2_PC(i,:,:)= data.GT_participation_coeff;
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group2_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group2_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group2_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group2_GE_rand(i,:,:,:)= data_rand.GT_global_eff;
    Group2_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group2_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
end
%%
Group2_CC_50=Group2_CC(:,1:20,:);                          
Group2_CC_rand_squ = squeeze(mean(Group2_CC_rand,3));        
AvgGroup2_CC=mean(mean(Group2_CC_50,3));

sparsity_CC_Group2_50 = (mean(Group2_CC_50,3));
sparsity_CC_rand_Group2 = mean(Group2_CC_rand_squ,3);
sparsity_CC_rand_Group2_50 = sparsity_CC_rand_Group2 (:,1:20);

Group2_PL_50=Group2_PL(:,:,1:20);
AvgGroup2_PL=squeeze(mean(Group2_PL_50));

Sparsity_PL_Group2_50=squeeze(Group2_PL_50);
Sparsity_PL_Group2_rand_squ = squeeze(mean(Group2_PL_rand,3)); 
Sparsity_PL_Group2_rand_50 = Sparsity_PL_Group2_rand_squ (:,1:20);
Group2_Degree_50=Group2_Degree(:,1:20,:);                             
AvgGroup2_Degree=mean(mean(Group2_Degree_50,3));
sparsity_Degree_Group2_50 = (mean(Group2_Degree_50,3));

Group2_PC_50=Group2_PC(:,1:20,:);                          
Group2_PC_rand_squ = squeeze(mean(Group2_PC_rand,3));        
AvgGroup2_PC=mean(mean(Group2_PC_50,3));

sparsity_PC_Group2_50 = (mean(Group2_PC_50,3));
sparsity_PC_rand_Group2 = mean(Group2_PC_rand_squ,3);
sparsity_PC_rand_Group2_50 = sparsity_PC_rand_Group2 (:,1:20);
%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            sparsity_CC_normalised_Group2(i,j) = sparsity_CC_Group2_50(i,j)/sparsity_CC_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            sparsity_PL_normalised_Group2(i,j) = Sparsity_PL_Group2_50(i,j)/Sparsity_PL_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            SmallWorldNess_Group2(i,j) = sparsity_CC_normalised_Group2(i,j)/sparsity_PL_normalised_Group2(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            sparsity_PC_normalised_Group2(i,j) = sparsity_PC_Group2_50(i,j)/sparsity_PC_rand_Group2_50(i,j); 
    end
end

%%
Group2_LE_50=Group2_LE(:,1:20,:);                          
Group2_LE_rand_squ = squeeze(mean(Group2_LE_rand,3));        
AvgGroup2_LE=mean(mean(Group2_LE_50,3));

sparsity_LE_Group2_50 = (mean(Group2_LE_50,3));
sparsity_LE_rand_Group2 = mean(Group2_LE_rand_squ,3);
sparsity_LE_rand_Group2_50 = sparsity_LE_rand_Group2 (:,1:20);

Group2_GE_50=Group2_GE(:,1:20);
AvgGroup2_GE=squeeze(mean(Group2_GE_50));

Sparsity_GE_Group2_50=squeeze(Group2_GE_50);
Sparsity_GE_Group2_rand_squ = squeeze(mean(Group2_GE_rand,3)); 
Sparsity_GE_Group2_rand_50 = Sparsity_GE_Group2_rand_squ (:,1:20); 

%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            sparsity_LE_normalised_Group2(i,j) = sparsity_LE_Group2_50(i,j)/sparsity_LE_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:20; 
            sparsity_GE_normalised_Group2(i,j) = Sparsity_GE_Group2_50(i,j)/Sparsity_GE_Group2_rand_50(i,j); 
    end
end

%% Absolute CC and PC ploting
figure (1)
y1_CC = mean (sparsity_CC_Group1_50);
z1_CC = std (sparsity_CC_Group1_50)/sqrt (length (sparsity_CC_Group1_50)); 
errorbar (y1_CC,z1_CC, 'b'); grid on; 
hold on
y2_CC = mean (sparsity_CC_Group2_50);
z2_CC = std (sparsity_CC_Group2_50)/sqrt (length (sparsity_CC_Group2_50));
errorbar (y2_CC,z2_CC, 'r'); grid on; 
title('Absolute Clustering Coeficent')

figure (2)
y1_PC = mean (sparsity_PC_Group1_50);
z1_PC = std (sparsity_PC_Group1_50)/sqrt (length (sparsity_PC_Group1_50)); 
errorbar (y1_PC,z1_PC, 'b'); grid on; 
hold on
y2_PC = mean (sparsity_PC_Group2_50);
z2_PC = std (sparsity_PC_Group2_50)/sqrt (length (sparsity_PC_Group2_50));
errorbar (y2_PC,z2_PC, 'r'); grid on; 
title('Absolute Participation Coeficent')
%% Random CC and PC ploting
figure (3)
y1_CC_rand = mean (sparsity_CC_rand_Group1_50);
z1_CC_rand = std (sparsity_CC_rand_Group1_50)/sqrt (length (sparsity_CC_rand_Group1_50)); 
errorbar (y1_CC_rand,z1_CC_rand, 'b'); grid on; 
hold on
y2_CC_rand = mean (sparsity_CC_rand_Group2_50);
z2_CC_rand = std (sparsity_CC_rand_Group2_50)/sqrt (length (sparsity_CC_rand_Group2_50));
errorbar (y2_CC_rand,z2_CC_rand, 'r'); grid on; 
title('Random Clustering Coeficent')

figure (4)
y1_PC_rand = mean (sparsity_PC_rand_Group1_50);
z1_PC_rand = std (sparsity_PC_rand_Group1_50)/sqrt (length (sparsity_PC_rand_Group1_50)); 
errorbar (y1_PC_rand,z1_PC_rand, 'b'); grid on; 
hold on
y2_PC_rand = mean (sparsity_PC_rand_Group2_50);
z2_PC_rand = std (sparsity_PC_rand_Group2_50)/sqrt (length (sparsity_PC_rand_Group2_50));
errorbar (y2_PC_rand,z2_PC_rand, 'r'); grid on; 
title('Random Participation Coeficent')
%%
figure (13), plot(sparsity_CC_normalised_Group1')
figure (14), plot(sparsity_CC_normalised_Group2')
figure (15), plot(sparsity_PC_normalised_Group1')
figure (16), plot(sparsity_PC_normalised_Group2')
%% %Ploting Normalized CC, PC, PL, GE, LE Images

figure (5)
y1_CC = mean (sparsity_CC_normalised_Group1);
z1_CC = std (sparsity_CC_normalised_Group1)/sqrt (length (sparsity_CC_normalised_Group1)); 
errorbar (y1_CC,z1_CC, 'b'); grid on; 
hold on
y2_CC = mean (sparsity_CC_normalised_Group2);
z2_CC = std (sparsity_CC_normalised_Group2)/sqrt (length (sparsity_CC_normalised_Group2));
errorbar (y2_CC,z2_CC, 'r'); grid on; 
title('Normalised Clustering Coeficent')

figure (6)
hold on
y1_PL = mean (sparsity_PL_normalised_Group1);
z1_PL = std (sparsity_PL_normalised_Group1)/sqrt (length (sparsity_PL_normalised_Group1)); 
errorbar (y1_PL,z1_PL, 'b'); grid on; 
hold on
y2_PL = mean (sparsity_PL_normalised_Group2);
z2_PL = std (sparsity_PL_normalised_Group2)/sqrt (length (sparsity_PL_normalised_Group2));
errorbar (y2_PL,z2_PL, 'r'); grid on;
title('Normalised Path Length')

figure (7)
y1_SW = mean (SmallWorldNess_Group1);
z1_SW = std (SmallWorldNess_Group1)/sqrt (length (SmallWorldNess_Group1));
errorbar (y1_SW,z1_SW, 'b'); grid on;
hold on
y2_SW = mean (SmallWorldNess_Group2);
z2_SW = std (SmallWorldNess_Group2)/sqrt (length (SmallWorldNess_Group2));
errorbar (y2_SW,z2_SW, 'r'); grid on;
title('Small Worldness')

figure (8)
y1_CC = mean (sparsity_LE_normalised_Group1);
z1_CC = std (sparsity_LE_normalised_Group1)/sqrt (length (sparsity_LE_normalised_Group1)); 
errorbar (y1_CC,z1_CC, 'b'); grid on; 
hold on
y2_CC = mean (sparsity_LE_normalised_Group2);
z2_CC = std (sparsity_LE_normalised_Group2)/sqrt (length (sparsity_LE_normalised_Group2));
errorbar (y2_CC,z2_CC, 'r'); grid on; 
title('Normalised Local Eficency')

figure (9)
hold on
y1_PL = mean (sparsity_GE_normalised_Group1);
z1_PL = std (sparsity_GE_normalised_Group1)/sqrt (length (sparsity_GE_normalised_Group1)); 
errorbar (y1_PL,z1_PL, 'b'); grid on; 
hold on
y2_PL = mean (sparsity_GE_normalised_Group2);
z2_PL = std (sparsity_GE_normalised_Group2)/sqrt (length (sparsity_GE_normalised_Group2));
errorbar (y2_PL,z2_PL, 'r'); grid on;
title('Normalised Global Eficency')
%
figure (10)
y1_PC = mean (sparsity_PC_normalised_Group1);
z1_PC = std (sparsity_PC_normalised_Group1)/sqrt (length (sparsity_PC_normalised_Group1)); 
errorbar (y1_PC,z1_PC, 'b'); grid on; 
hold on
y2_PC = mean (sparsity_PC_normalised_Group2);
z2_PC = std (sparsity_PC_normalised_Group2)/sqrt (length (sparsity_PC_normalised_Group2));
errorbar (y2_PC,z2_PC, 'r'); grid on; 
title('Normalised Participation Coeficent')

%% --------------t-stats between Group2 and Group1 sparsity level------------%

for i = 1:20
    [h_CC1(i),p_CC1(i)] = ttest2(sparsity_CC_normalised_Group2(:,i),sparsity_CC_normalised_Group1(:,i),0.05,'left');
end

h_CC1
p_CC1

for i = 1:20
    [h_PL1(i),p_PL1(i)] = ttest2(sparsity_PL_normalised_Group2(:,i),sparsity_PL_normalised_Group1(:,i),0.05,'right');
end
h_PL1
p_PL1
%%
for i = 1:20
    [h_SW1(i),p_SW1(i)] = ttest2(SmallWorldNess_Group2(:,i),SmallWorldNess_Group1(:,i),0.05,'left');
end

h_SW1
p_SW1
%%
for i = 1:20
    [h_PC1(i),p_PC1(i)] = ttest2(sparsity_PC_normalised_Group2(:,i),sparsity_PC_normalised_Group1(:,i),0.05,'left');
end

h_PC1
p_PC1
%%
%----------------------%% Brain resion significant Computations for CC %-------------%


sparsity_CC_Group1_ROI = squeeze(mean(Group1_CC_50,2)); 
sparsity_CC_Group2_ROI = squeeze(mean(Group2_CC_50,2));

sparsity_CC_rand_Group1_ROI = mean(Group1_CC_rand_squ,2);
sparsity_CC_rand_Group2_ROI = mean(Group2_CC_rand_squ,2);



for i = 1:length(SUBJlist_Group1); 
    for j = 1:size(sparsity_PC_Group1_ROI,2); 
            sparsity_CC_normalised_Group1_ROI(i,j) = sparsity_CC_Group1_ROI(i,j)/sparsity_CC_rand_Group1_ROI(i,j);
    end
end



for i = 1:length(SUBJlist_Group2); 
    for j = 1:size(sparsity_PC_Group2_ROI,2); 
            sparsity_CC_normalised_Group2_ROI(i,j) = sparsity_CC_Group2_ROI(i,j)/sparsity_CC_rand_Group2_ROI(i,j); 
    end
end

%%
%----------------------%% Brain resion significant Computations for PC %-------------%


sparsity_PC_Group1_ROI = squeeze(mean(Group1_PC_50,2)); 
sparsity_PC_Group2_ROI = squeeze(mean(Group2_PC_50,2));

sparsity_PC_rand_Group1_ROI = mean(Group1_PC_rand_squ,2);
sparsity_PC_rand_Group2_ROI = mean(Group2_PC_rand_squ,2);



for i = 1:length(SUBJlist_Group1); 
    for j = 1:size(sparsity_PC_Group1_ROI,2); 
            sparsity_PC_normalised_Group1_ROI(i,j) = sparsity_PC_Group1_ROI(i,j)/sparsity_PC_rand_Group1_ROI(i,j);
    end
end



for i = 1:length(SUBJlist_Group2); 
    for j = 1:size(sparsity_PC_Group2_ROI,2); 
            sparsity_PC_normalised_Group2_ROI(i,j) = sparsity_PC_Group2_ROI(i,j)/sparsity_PC_rand_Group2_ROI(i,j); 
    end
end


%% Brain resion significant Computations for PC

for i = 1:31
    [h_CC_ROI_Normalised(i),p_CC_ROI_Normalised(i)] = ttest2(sparsity_CC_normalised_Group1_ROI(:,i),sparsity_CC_normalised_Group2_ROI(:,i),0.02,'left');
end

h_CC_ROI_Normalised
p_CC_ROI_Normalised
%%
for i = 1:31
    [h_PC_ROI_Normalised(i),p_PC_ROI_Normalised(i)] = ttest2(sparsity_PC_normalised_Group1_ROI(:,i),sparsity_PC_normalised_Group2_ROI(:,i),0.02,'right');
end

h_PC_ROI_Normalised   
p_PC_ROI_Normalised
%%[p,h]=fdr(p_CC_ROI_Normalised,0.05);
%p

%%  Graphical Plot od GT measures %%%
tpz_cc_control = trapz(sparsity_CC_normalised_Group1(:,1:12),2)/12;
tpz_cc_Patient = trapz(sparsity_CC_normalised_Group2(:,1:12),2)/12;
tpz_pc_control = trapz(sparsity_PC_normalised_Group1(:,1:20),2)/20;
tpz_pc_Patient = trapz(sparsity_PC_normalised_Group2(:,1:20),2)/20;
tpz_pl_control = trapz(sparsity_PL_normalised_Group1(:,1:12),2)/12;
tpz_pl_Patient = trapz(sparsity_PL_normalised_Group2(:,1:12),2)/12;
tpz_sw_control = trapz(SmallWorldNess_Group1(:,1:12),2)/12;
tpz_sw_Patient = trapz(SmallWorldNess_Group2(:,1:12),2)/12;
tpz_le_control = trapz(sparsity_LE_normalised_Group1(:,1:12),2)/12;
tpz_le_Patient = trapz(sparsity_LE_normalised_Group2(:,1:12),2)/12;
tpz_ge_control = trapz(sparsity_GE_normalised_Group1(:,1:12),2)/12;
tpz_ge_Patient = trapz(sparsity_GE_normalised_Group2(:,1:12),2)/12;

%tpz_GT = [tpz_cc_control tpz_cc_Patient tpz_pc_control tpz_pc_Patient tpz_pl_control tpz_pl_Patient tpz_sw_control tpz_sw_Patient tpz_le_control tpz_le_Patient tpz_ge_control tpz_ge_Patient]
tpz_GT_control = [tpz_cc_control tpz_pc_control tpz_pl_control tpz_sw_control tpz_le_control tpz_ge_control];
tpz_GT_Patient = [tpz_cc_Patient tpz_pc_Patient tpz_pl_Patient tpz_sw_Patient tpz_le_Patient tpz_ge_Patient];

% % % [p,h]=ttest2(tpz_GT(:,1),tpz_GT(:,2), 0.05,'left')
% % % 
% % % 


%%
%% Network Wise Normalised PC and CC calculation
SN=1; FP=4; DMN=3; SubCor=4; SMN=5; Visual1=6; Visual2=7; Visual3=8; Cer=9;
Grp_1 = ones(length(SUBJlist_Group1),1);
Grp_2 = ones(length(SUBJlist_Group2),1)*2;
Grp = [Grp_1; Grp_2];
%SN=2; FP=4; DMN=3; SubCor=1; SMN=6; Visual1=5; % Docenbach % load('J:\EEG\FrontoParietalDOC\Processes\ZGT\DocenbachNetworkNo.mat')
load('J:\fMRI\FrontoParietalDOC\GT_DFC\ShnenNetworkNo.mat')
Grp_PC_SN = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SN)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SN)),2)];
figure(115); notBoxPlot(Grp_PC_SN,Grp,0.5,'patch',ones(length(Grp_PC_SN),1));
Grp_PC_FP = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,FP)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,FP)),2); mean(sparsity_PC_normalised_Group3_ROI(:, ismember(ShnenNetworkNo,FP)),2)];
figure(116); notBoxPlot(Grp_PC_FP,Grp,0.5,'patch',ones(length(Grp_PC_FP),1));
Grp_PC_DMN = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,DMN)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,DMN)),2); mean(sparsity_PC_normalised_Group3_ROI(:, ismember(ShnenNetworkNo,DMN)),2)];
figure(117); notBoxPlot(Grp_PC_DMN,Grp,0.5,'patch',ones(length(Grp_PC_DMN),1));
Grp_PC_SubCor = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SubCor)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SubCor)),2); mean(sparsity_PC_normalised_Group3_ROI(:, ismember(ShnenNetworkNo,SubCor)),2)];
figure(118); notBoxPlot(Grp_PC_SubCor,Grp,0.5,'patch',ones(length(Grp_PC_SubCor),1));
Grp_PC_SMN = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SMN)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SMN)),2); mean(sparsity_PC_normalised_Group3_ROI(:, ismember(ShnenNetworkNo,SMN)),2)];
figure(119); notBoxPlot(Grp_PC_SMN,Grp,0.5,'patch',ones(length(Grp_PC_SMN),1));
Grp_PC_Visual1 = [mean(sparsity_PC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,Visual1)),2); mean(sparsity_PC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,Visual1)),2); mean(sparsity_PC_normalised_Group3_ROI(:, ismember(ShnenNetworkNo,Visual1)),2)];
figure(120); notBoxPlot(Grp_PC_Visual1,Grp,0.5,'patch',ones(length(Grp_PC_Visual1),1));
%% ---
SN=1; FP=4; DMN=3; SubCor=4; SMN=5; Visual1=6; Visual2=7; Visual3=8; Cer=9;
Grp_1 = ones(length(SUBJlist_Group1),1);
Grp_2 = ones(length(SUBJlist_Group2),1)*2;
Grp = [Grp_1; Grp_2];
%SN/Auditory=2; FP=4; DMN=3; SubCor=1; SMN=6; Visual1=5; % Docenbach % load('J:\EEG\FrontoParietalDOC\Processes\ZGT\DocenbachNetworkNo.mat')
load('J:\fMRI\FrontoParietalDOC\GT_DFC\ShnenNetworkNo.mat')
Grp_CC_SN = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SN)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SN)),2)];
figure(125); notBoxPlot(Grp_CC_SN,Grp,0.5,'patch',ones(length(Grp_CC_SN),1));
title('Auditory Network');
Grp_CC_FP = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,FP)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,FP)),2)];
figure(126); notBoxPlot(Grp_CC_FP,Grp,0.5,'patch',ones(length(Grp_CC_FP),1));
title('Fronto-Parietalt Network');
Grp_CC_DMN = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,DMN)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,DMN)),2)];
figure(127); notBoxPlot(Grp_CC_DMN,Grp,0.5,'patch',ones(length(Grp_CC_DMN),1));
title('DMN Network');
Grp_CC_SubCor = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SubCor)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SubCor)),2)];
figure(128); notBoxPlot(Grp_CC_SubCor,Grp,0.5,'patch',ones(length(Grp_CC_SubCor),1));
title('Sub Cortical Network');
Grp_CC_SMN = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,SMN)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,SMN)),2)];
figure(129); notBoxPlot(Grp_CC_SMN,Grp,0.5,'patch',ones(length(Grp_CC_SMN),1));
title('Sensory Motor Network');
Grp_CC_Visual1 = [mean(sparsity_CC_normalised_Group1_ROI(:, ismember(ShnenNetworkNo,Visual1)),2); mean(sparsity_CC_normalised_Group2_ROI(:, ismember(ShnenNetworkNo,Visual1)),2)];
figure(130); notBoxPlot(Grp_CC_Visual1,Grp,0.5,'patch',ones(length(Grp_CC_Visual1),1));
title('Visual Network');