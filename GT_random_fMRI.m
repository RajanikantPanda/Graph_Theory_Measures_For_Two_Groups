%  Author: Dr. Rajanikant Panda (Email:bk.bme.rajanikant@gmail.com).
%  Description: Computes graph theory measures for random matrix, including:
%  degree, clustering coefficient, participation coefficient, path length, small-worldness, efficiency, and modularity.clc
%  this can be used to compute normalised matrices, that is Absolute/Random

clc
clear
path='J:\fMRI\TINNITUS\Patients\';
cd(path)
SUBJlist=dir('dataanalysis_Pnt*');
SUBJlist=SUBJlist(1:length(SUBJlist));
%%
for i=1:length(SUBJlist)
    SUBJname=SUBJlist(i).name;
    path1=([path SUBJname])
    cd(path1);
    %% 
    filelist= ([ SUBJname ]);  
        
    for i1=1
        filename=filelist; 
        prefix_name=filename(14:end);
                
        final_data=load(filename);
        final_data=final_data.y_roi_regressed_filtered;

        GT_corr_data=corr(final_data);
        GT_corr_data_abs = GT_corr_data;
        GT_corr_data_abs(GT_corr_data_abs < 0) = 0;
        chanlocs = size(GT_corr_data,1); % No  of ROIs                       
        %% Thresholding
        sparsity_val=0.01:0.025:0.5;
        for i2=1:20
            %% %%%Network properties/ network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
                       
            for random_number=1:50
                             
                random_network=randmio_und(corr_data_thr,5);
                %%
                GT_corr_data_rand_thr(i2,random_number,:,:)=random_network;% asign to different array
                
                GT_degree_rand(i2,random_number,:)=degrees_und(random_network); % calculate degree
                
                GT_clust_coeff_rand(i2,random_number,:)=clustering_coef_bu(random_network);%% calculate clustering coeff
                                
                GT_local_eff_rand(i2,random_number,:)=efficiency_bin(random_network,1); % local efficiency
                
                GT_global_eff(i2,random_number,:)=efficiency_bin(random_network,0); % global efficiency
                
                GT_distance_matrix_rand(i2,random_number,:,:)=distance_bin(random_network); % distance matrix
                
                GT_path_length_rand(i2,random_number) =charpath(squeeze(GT_distance_matrix_rand(i2,random_number,:,:)),1,0); % path length
            
            %%%%% Participation coefficient and modspan
                param.heuristic=50;
                for i = 1:param.heuristic
                    [Ci, allQ(i2,random_number,i)] = community_louvain(random_network);
                     allCi(i2,random_number,i,:) = Ci;
                     allpc(i2,random_number,i,:) = participation_coef(random_network,Ci); 
                end
                GT_modularity_rand(i2,random_number)= mean(allQ(i2,random_number,:));  % modularity
                GT_community_structure_rand(i2,random_number,1:chanlocs) = squeeze(allCi(i2,random_number,1,:)); % community structure
                GT_participation_coeff_rand(i2,random_number,1:chanlocs) = mean(squeeze(allpc(i2,random_number,:,:)));  %participation coefficient
                %%%%%    
            end
        end
        varname=([SUBJname '_RAND'])
        save(varname);
    end
    cd ..
end

