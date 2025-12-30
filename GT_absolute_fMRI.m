%  Author: Dr. Rajanikant Panda (Email:bk.bme.rajanikant@gmail.com) 
%  Description: Computes absolute graph theory measures, including:
%  degree, clustering coefficient, participation coefficient, path length, small-worldness, efficiency, and modularity.clc

clear
path='J:\fMRI\TINNITUS\Control\';
cd(path)
SUBJlist=dir('dataanalysis_Cnt*');
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
        chanlocs = size(final_data,2); % No  of ROIs      
        %% Thresholding
        sparsity_val=0.01:0.025:0.5; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)    
            %% %%%Network properties/ network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
            
            %%
            GT_corr_data_thr(i2,:,:)=corr_data_thr;% asign to different array
            
            GT_degree(i2,:)=degrees_und(corr_data_thr); % calculate degree
            
            GT_clust_coeff(i2,:)=clustering_coef_bu(corr_data_thr);%% calculate clustering coeff
            
            GT_local_eff(i2,:)=efficiency_bin(corr_data_thr,1); % local efficiency
            
            GT_global_eff(i2,:)=efficiency_bin(corr_data_thr,0); % global efficiency
            
            GT_distance_matrix(i2,:,:)=distance_bin(corr_data_thr); % distance matrix
            
            GT_path_length(i2)=charpath(squeeze(GT_distance_matrix(i2,:,:)),1,0); % path length
            %path_length(s)  = charpath(distance_matrix(s,:,:),1,0); % path length
              
            %%%%% Participation coefficient and modspan
            param.heuristic=50;
            
            for i = 1:param.heuristic
                [Ci, allQ(i2,i)] = community_louvain(corr_data_thr);
            
                allCi(i2,i,:) = Ci;
   
                allpc(i2,i,:) = participation_coef(corr_data_thr,Ci); 
            end
        
            GT_modularity(i2)= mean(allQ(i2,:));  % modularity
            GT_community_structure(i2,1:chanlocs) = squeeze(allCi(i2,1,:)); % community structure
            GT_participation_coeff(i2,1:chanlocs) = mean(squeeze(allpc(i2,:,:)));  % participation coefficient
            
            %%%%%
            
        end
        
        varname=([SUBJname '_ABS'])
        save(varname);
    end
    cd ..

end
