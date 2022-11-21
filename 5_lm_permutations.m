% Perform permutation analyses to obtain null distribution and non-parametric p values for linear  models (lm) associating FC similarity and behavioural measure changes

% This script assumes the lm models have been run and saved for FC similarity and behavioural change associations to use as inputs
% This script also assumes the lme models have been run and saved for FC similarity and behavioural scores to use as inputs

% All parts that need to be updated with specific information are indicated with TO_BE_SPECIFIED

clear;clc % clear workspace variables and command window

n_perm = TO_BE_SPECIFIED; % number of permutations to run

load('TO_BE_SPECIFIED.mat') % input file containing FC similarity values
load('TO_BE_SPECIFIED.mat') % input file containing lme models
clear all_behav

network = {'TO_BE_SPECIFIED'}; % network names
behav = {'TO_BE_SPECIFIED'}; % behavioural scores used

network_lme_tbl = inter_lme_mdl{length(network)}.Variables; 
behav_lme_tbl = behav_lme_mdl{length(behav)}.Variables;
lme_tbl = [behav_lme_tbl network_lme_tbl(:,6:end)]; % store behaviour and FC similarity scores in same table (without repeating demographic variables)

grp_perm_idx = NaN(length(behav_lme_tbl.group)/2,n_perm); % to store shuffled group labels
time_perm = randi(2,length(behav_lme_tbl.group)/2,n_perm); % flip timepoint labels or not (50% chance) - for pre
time_perm = categorical([time_perm;3-time_perm]); % update post accordingly

for perm_num = 1:n_perm
    grp_perm_idx(:,perm_num) = randperm(length(lme_tbl.group)/2); % shuffle group labels for each iteration of permutations
	
	% update variables table with permuted group and timepoint labels 
    lme_tbl_perm = lme_tbl;
    lme_tbl_perm.group = [lme_tbl.group(grp_perm_idx(:,perm_num));lme_tbl.group(grp_perm_idx(:,perm_num))];
    lme_tbl_perm.timepoint = time_perm(:,perm_num);
    lme_tbl_perm = sortrows(lme_tbl_perm,'timepoint');
	
	% collect demographics to use in lm
	lm_tbl_perm = lme_tbl_perm(1:size(lme_tbl,1)/2,[1 3:5]); % excluding timepoint (not used in lm)
    
    for behav_num = 1:length(behav) % get permuted behavioural changes
        behav_data_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.Properties.VariableNames,behav{behav_num}))); % get column idx for each behaviour
        for subj_num = 1:size(lme_tbl,1)/2
            subj_data_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.subject,subj_name{subj_num}))); % get idx for the subject (two rows per subject)
            behav_diff(subj_num,behav_num) = lme_tbl_perm{subj_data_idx(2),behav_data_idx} - lme_tbl_perm{subj_data_idx(1),behav_data_idx};  % behavioural change
        end
        lm_tbl_perm = [lm_tbl_perm array2table(behav_diff(:,behav_num),'VariableNames',{[behav{behav_num} '_diff']})]; % add behavioural change to lm table
    end
	
    for network_num = 1:length(network) % get permuted FC similarity changes
        network_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.Properties.VariableNames,[network{network_num} '_intra_fc_sim']))); % get column idx for intranework FC similarity
        for subj_num = 1:size(lme_tbl,1)/2
            subj_data_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.subject,subj_name{subj_num}))); % get idx for the subject (two rows per subject)
            fcsim_diff(subj_num,behav_num) = lme_tbl_perm{subj_data_idx(2),network_idx} - lme_tbl_perm{subj_data_idx(1),network_idx}; % FC similarity change
        end
        lm_tbl_perm = [lm_tbl_perm array2table(fcsim_diff(:,behav_num),'VariableNames',{[network{network_num} '_intra_fc_sim_diff']})]; % add intranetwork FC similarity change to lm table
		
        network_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.Properties.VariableNames,[network{network_num} '_inter_fc_sim']))); % get column idx for internework FC similarity
        for subj_num = 1:size(lme_tbl,1)/2
            subj_data_idx = find(~cellfun('isempty',strfind(lme_tbl_perm.subject,subj_name{subj_num}))); % get idx for the subject (two rows per subject)
            fcsim_diff(subj_num,behav_num) = lme_tbl_perm{subj_data_idx(2),network_idx} - lme_tbl_perm{subj_data_idx(1),network_idx}; % FC similarity change
        end
        lm_tbl_perm = [lm_tbl_perm array2table(fcsim_diff(:,behav_num),'VariableNames',{[network{network_num} '_inter_fc_sim_diff']})]; % add internetwork FC similarity change to lm table
    end
    
	% do permutation lm to get null distribution
    for behav_num = 1:length(behav)
        for network_num = 1:length(network)
            
			% permutation lm with group term
            intra_lm_mdl_perm{network_num,behav_num,perm_num} = fitlm(lm_tbl_perm,[behav{behav_num} '_diff~age+gender+group*' network{network_num} '_intra_fc_sim_diff']);
            inter_lm_mdl_perm{network_num,behav_num,perm_num} = fitlm(lm_tbl_perm,[behav{behav_num} '_diff~age+gender+group*' network{network_num} '_inter_fc_sim_diff']);
            
			% collect main effect values from lm
			idx = find(strcmp(intra_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,[network{network_num} '_intra_fc_sim_diff'])); % association
            intra_FCsim_beta_null(network_num,behav_num,perm_num) = intra_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
            idx = find(strcmp(intra_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,['group_2:' network{network_num} '_intra_fc_sim_diff'])); % interaction with group
            intra_int_beta_null(network_num,behav_num,perm_num) = intra_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
            idx = find(strcmp(inter_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,[network{network_num} '_inter_fc_sim_diff'])); % association
            inter_FCsim_beta_null(network_num,behav_num,perm_num) = inter_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
            idx = find(strcmp(inter_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,['group_2:' network{network_num} '_inter_fc_sim_diff'])); % interaction with group
            inter_int_beta_null(network_num,behav_num,perm_num) = inter_lm_mdl_perm{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
            
			% permutation lm without group term			
			intra_lm_mdl_perm_c{network_num,behav_num,perm_num} = fitlm(lm_tbl_perm,[behav{behav_num} '_diff~age+gender+' network{network_num} '_intra_fc_sim_diff']);
            inter_lm_mdl_perm_c{network_num,behav_num,perm_num} = fitlm(lm_tbl_perm,[behav{behav_num} '_diff~age+gender+' network{network_num} '_inter_fc_sim_diff']);
                       
			% collect main effect values from lm
            idx = find(strcmp(intra_lm_mdl_perm_c{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,[network{network_num} '_intra_fc_sim_diff'])); % association
            intra_FCsim_beta_null_c(network_num,behav_num,perm_num) = intra_lm_mdl_perm_c{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
            idx = find(strcmp(inter_lm_mdl_perm_c{network_num,behav_num,perm_num}.Coefficients.Properties.RowNames,[network{network_num} '_inter_fc_sim_diff'])); % association
            inter_FCsim_beta_null_c(network_num,behav_num,perm_num) = inter_lm_mdl_perm_c{network_num,behav_num,perm_num}.Coefficients.Estimate(idx);
        end
    end
    clear fcsim_diff behav_diff
end

%% get p values based on null distributions from permutation lm

% headers for output (p values) to be collected - lm with group term
network_perm_output = {};
network_perm_output(1) = {'behav,network,network_type,FCsim_p,interaction_p'};

load('TO_BE_SPECIFIED.mat') % input file containing lm (with group term) associating FC similarity change with behavioural change

j = 2;
for behav_num = 1:length(behav)
    for network_num = 1:length(network)
	
		% [intranetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
        idx = find(strcmp(intra_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,[network{network_num} '_intra_fc_sim_diff'])); % association
        if intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            time_p = sum(intra_FCsim_beta_null(network_num,behav_num,:)<=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            time_p = sum(intra_FCsim_beta_null(network_num,behav_num,:)>=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        idx = find(strcmp(intra_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,['group_2:' network{network_num} '_intra_fc_sim_diff'])); % interaction with group
        if intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            interaction_p = sum(intra_int_beta_null(network_num,behav_num,:)<=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            interaction_p = sum(intra_int_beta_null(network_num,behav_num,:)>=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        network_perm_output(j) = strcat(behav(behav_num),',',network(network_num),',intra,',num2str(time_p),',',num2str(interaction_p)); % store permutation p values
        j = j+1;

		% [internetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
        idx = find(strcmp(inter_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,[network{network_num} '_inter_fc_sim_diff'])); % association
        if inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            time_p = sum(inter_FCsim_beta_null(network_num,behav_num,:)<=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            time_p = sum(inter_FCsim_beta_null(network_num,behav_num,:)>=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        idx = find(strcmp(inter_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,['group_2:' network{network_num} '_inter_fc_sim_diff'])); % interaction with group
        if inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            interaction_p = sum(inter_int_beta_null(network_num,behav_num,:)<=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            interaction_p = sum(inter_int_beta_null(network_num,behav_num,:)>=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        network_perm_output(j) = strcat(behav(behav_num),',',network(network_num),',inter,',num2str(time_p),',',num2str(interaction_p)); % store permutation p values
        j = j+1;
    end
end

% headers for output (p values) to be collected - lm without group term
network_perm_output_c = {};
network_perm_output_c(1) = {'behav,network,network_type,FCsim_p'};

load('TO_BE_SPECIFIED.mat') % input file containing lm (no group term) associating FC similarity change with behavioural change

j = 2;
for behav_num = 1:length(behav)
    for network_num = 1:length(network)
	
		% [intranetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
        idx = find(strcmp(intra_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,[network{network_num} '_intra_fc_sim_diff'])); % association
        if intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            time_p = sum(intra_FCsim_beta_null_c(network_num,behav_num,:)<=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            time_p = sum(intra_FCsim_beta_null_c(network_num,behav_num,:)>=intra_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        network_perm_output_c(j) = strcat(behav(behav_num),',',network(network_num),',intra,',num2str(time_p)); % store permutation p values
        j=j+1;
		
		% [internetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
        idx = find(strcmp(inter_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,[network{network_num} '_inter_fc_sim_diff'])); % association
        if inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) < 0
            time_p = sum(inter_FCsim_beta_null_c(network_num,behav_num,:)<=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        elseif inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx) > 0
            time_p = sum(inter_FCsim_beta_null_c(network_num,behav_num,:)>=inter_lm_mdl{network_num,behav_num}.Coefficients.Estimate(idx))/n_perm;
        end
        network_perm_output_c(j) = strcat(behav(behav_num),',',network(network_num),',inter,',num2str(time_p)); % store permutation p values        
        j = j+1;
    end
end

% readable outputs saved as csv files
writetable(cell2table(network_perm_output'),'TO_BE_SPECIFIED.csv'); % with group term
writetable(cell2table(network_perm_output_c'),'TO_BE_SPECIFIED.csv'); % without group term

% null distributions from permutations and permutation indices saved in mat file
save('TO_BE_SPECIFIED.mat','grp_perm_idx','time_perm','intra_FCsim_beta_null','intra_int_beta_null','inter_FCsim_beta_null','inter_int_beta_null'); % with group term
save('TO_BE_SPECIFIED.mat','grp_perm_idx','time_perm','intra_FCsim_beta_null_c','inter_FCsim_beta_null_c'); % without group term
