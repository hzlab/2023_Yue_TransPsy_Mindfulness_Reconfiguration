% Perform permutation analyses to obtain null distribution and non-parametric p values for linear mixed models (lme) on FC similarity and behavioural measures

% This script assumes the lme models have been run and saved for FC similarity and behavioural scores to use as inputs
% All parts that need to be updated with specific information are indicated with TO_BE_SPECIFIED

clear;clc % clear workspace variables and command window

n_perm = TO_BE_SPECIFIED; % number of permutations to run

load('TO_BE_SPECIFIED.mat') % file containing lme for FC similarity and behaviour
network = {'TO_BE_SPECIFIED'}; % network names
behav = {'TO_BE_SPECIFIED'}; % behavioural scores used

% get table of variables from loaded lme models
network_lme_tbl = inter_lme_mdl{length(network)}.Variables; 
behav_lme_tbl = behav_lme_mdl{length(behav)}.Variables;

grp_perm_idx = NaN(length(behav_lme_tbl.group)/2,n_perm); % to store shuffled group labels
time_perm = randi(2,length(behav_lme_tbl.group)/2,n_perm); % flip timepoint labels or not (50% chance) - for pre
time_perm = categorical([time_perm;3-time_perm]); % update post accordingly

for perm_num = 1:n_perm
    grp_perm_idx(:,perm_num) = randperm(length(network_lme_tbl.group)/2); % shuffle group labels for each iteration of permutations
	
    for network_num = 1:length(network) % do permutations for FC similarity
	
		% update variables table with permuted group and timepoint labels 
        network_lme_tbl_perm = network_lme_tbl;
        network_lme_tbl_perm.group = [network_lme_tbl.group(grp_perm_idx(:,perm_num));network_lme_tbl.group(grp_perm_idx(:,perm_num))];
        network_lme_tbl_perm.timepoint = time_perm(:,perm_num);

		% collect fixed effect values from models in each permutation
        intra_lme_mdl_perm{network_num,perm_num} = fitlme(network_lme_tbl_perm,[network{network_num} '_intra_fc_sim~baseline_age+gender+group*timepoint+(timepoint|subject)']); 
        idx = find(strcmp(intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'timepoint_2')); % main effect of time
        intra_time_beta_null(network_num,perm_num) = intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'group_2')); % main effect of group
        intra_group_beta_null(network_num,perm_num) = intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'timepoint_2:group_2')); % interaction effect
        intra_int_beta_null(network_num,perm_num) = intra_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx);
        inter_lme_mdl_perm{network_num,perm_num} = fitlme(network_lme_tbl_perm,[network{network_num} '_inter_fc_sim~baseline_age+gender+group*timepoint+(timepoint|subject)']); 
        idx = find(strcmp(inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'timepoint_2')); % main effect of time
        inter_time_beta_null(network_num,perm_num) = inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'group_2')); % main effect of group
        inter_group_beta_null(network_num,perm_num) = inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Name,'timepoint_2:group_2')); % interaction effect
        inter_int_beta_null(network_num,perm_num) = inter_lme_mdl_perm{network_num,perm_num}.Coefficients.Estimate(idx); 
        
        % collect predicted values and do pairwise comparisons (post-hoc) for each permutation 
        idx_pre = find(network_lme_tbl_perm.timepoint == '1');
        idx_post = find(network_lme_tbl_perm.timepoint == '2');
        idx_mind = find(network_lme_tbl_perm.group == '1'); % MBTI
        idx_slp = find(network_lme_tbl_perm.group == '2'); % SHEEP 
        pred_val = predict(intra_lme_mdl_perm{network_num,perm_num},network_lme_tbl_perm); % predicted values from permuted lme
        network_posthoc_null(network_num,1,1,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
        network_posthoc_null(network_num,1,2,perm_num) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
        network_posthoc_null(network_num,1,3,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
        network_posthoc_null(network_num,1,4,perm_num) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
        pred_val = predict(inter_lme_mdl_perm{network_num,perm_num},network_lme_tbl_perm); % predicted values from permuted lme
        network_posthoc_null(network_num,2,1,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
        network_posthoc_null(network_num,2,2,perm_num) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
        network_posthoc_null(network_num,2,3,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
        network_posthoc_null(network_num,2,4,perm_num) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
    end
    
    for behav_num = 1:length(behav) % do permutations for behavioural scores
		
		% update variables table with permuted group and timepoint labels 
        behav_lme_tbl_perm = behav_lme_tbl;
        behav_lme_tbl_perm.group = [behav_lme_tbl.group(grp_perm_idx(:,perm_num));behav_lme_tbl.group(grp_perm_idx(:,perm_num))];
        behav_lme_tbl_perm.timepoint = time_perm(:,perm_num);

		% collect fixed effect values from models in each permutation		
        behav_lme_mdl_perm{behav_num,perm_num} = fitlme(behav_lme_tbl_perm,[behav{behav_num} '~age+gender+group*timepoint+(timepoint|subject)']);
        idx = find(strcmp(behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Name,'timepoint_2')); % main effect of time 
        behav_time_beta_null(behav_num,perm_num) = behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Name,'group_2')); % main effect of group
        behav_group_beta_null(behav_num,perm_num) = behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Estimate(idx); 
        idx = find(strcmp(behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Name,'timepoint_2:group_2')); % interaction effect
        behav_int_beta_null(behav_num,perm_num) = behav_lme_mdl_perm{behav_num,perm_num}.Coefficients.Estimate(idx); 
        
        % collect predicted values and do pairwise comparisons (post-hoc) for each permutation 
        idx_pre = find(behav_lme_tbl_perm.timepoint == '1');
        idx_post = find(behav_lme_tbl_perm.timepoint == '2');
        idx_mind = find(behav_lme_tbl_perm.group == '1'); % MBTI
        idx_slp = find(behav_lme_tbl_perm.group == '2'); % SHEEP
		pred_val = predict(behav_lme_mdl_perm{behav_num,perm_num},behav_lme_tbl_perm); % predicted values from permuted lme
        behav_posthoc_null(behav_num,1,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
        behav_posthoc_null(behav_num,2,perm_num) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
        behav_posthoc_null(behav_num,3,perm_num) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
        behav_posthoc_null(behav_num,4,perm_num) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
    end
end
            
%% get p values based on null distributions from permutation lme

% headers for output (p values) to be collected
network_perm_output(1) = {'network,network_type,time_p,group_p,interaction_p'};
network_posthoc_output(1) = {'network,network_type,preM_postM,preS_postS,preM_preS,postM_postS'};
behav_perm_output(1) = {'behav,time_p,group_p,interaction_p'};
behav_posthoc_output(1) = {'behav,preM_postM,preS_postS,preM_preS,postM_postS'};

% collect p values for FC similarity
j = 2;
k = 2;
for network_num = 1:length(network) 

	% [intranetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
    idx = find(strcmp(intra_lme_mdl{network_num}.Coefficients.Name,'timepoint_2')); % for main effect of time
    if intra_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        time_p = sum(intra_time_beta_null(network_num,:)<=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif intra_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        time_p = sum(intra_time_beta_null(network_num,:)>=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(intra_lme_mdl{network_num}.Coefficients.Name,'group_2')); % for main effect of group
    if intra_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        group_p = sum(intra_group_beta_null(network_num,:)<=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif intra_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        group_p = sum(intra_group_beta_null(network_num,:)>=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(intra_lme_mdl{network_num}.Coefficients.Name,'timepoint_2:group_2')); % for interaction effect
    if intra_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        interaction_p = sum(intra_int_beta_null(network_num,:)<=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif intra_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        interaction_p = sum(intra_int_beta_null(network_num,:)>=intra_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    network_perm_output(j) = strcat(network(network_num),',intra,',num2str(time_p),',',num2str(group_p),',',num2str(interaction_p)); % store permutation p values
    j = j+1;
	
    % [internetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
    idx = find(strcmp(inter_lme_mdl{network_num}.Coefficients.Name,'timepoint_2')); % for main effect of time
    if inter_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        time_p = sum(inter_time_beta_null(network_num,:)<=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif inter_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        time_p = sum(inter_time_beta_null(network_num,:)>=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(inter_lme_mdl{network_num}.Coefficients.Name,'group_2')); % for main effect of group
    if inter_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        group_p = sum(inter_group_beta_null(network_num,:)<=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif inter_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        group_p = sum(inter_group_beta_null(network_num,:)>=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(inter_lme_mdl{network_num}.Coefficients.Name,'timepoint_2:group_2')); % for interaction effect
    if inter_lme_mdl{network_num}.Coefficients.Estimate(idx) < 0
        interaction_p = sum(inter_int_beta_null(network_num,:)<=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    elseif inter_lme_mdl{network_num}.Coefficients.Estimate(idx) > 0
        interaction_p = sum(inter_int_beta_null(network_num,:)>=inter_lme_mdl{network_num}.Coefficients.Estimate(idx))/n_perm;
    end
    network_perm_output(j) = strcat(network(network_num),',inter,',num2str(time_p),',',num2str(group_p),',',num2str(interaction_p)); % store permutation p values
    j = j+1;
    
    % [intranetwork] calculate observed posthoc (pairwise comparison) values
	idx_pre = find(network_lme_tbl.timepoint == '1');
    idx_post = find(network_lme_tbl.timepoint == '2');
    idx_mind = find(network_lme_tbl.group == '1'); % MBTI
    idx_slp = find(network_lme_tbl.group == '2'); % SHEEP
	pred_val = predict(intra_lme_mdl{network_num},network_lme_tbl);
    posthoc(1) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
    posthoc(2) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
    posthoc(3) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
    posthoc(4) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
	
	% [intranetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed posthoc values
    if posthoc(1) < 0
        posthoc_p1 = sum(network_posthoc_null(network_num,1,1,:)<=posthoc(1))/n_perm;
    elseif posthoc(1) > 0
        posthoc_p1 = sum(network_posthoc_null(network_num,1,1,:)>=posthoc(1))/n_perm;
    end
    if posthoc(2) < 0
        posthoc_p2 = sum(network_posthoc_null(network_num,1,2,:)<=posthoc(2))/n_perm;
    elseif posthoc(2) > 0
        posthoc_p2 = sum(network_posthoc_null(network_num,1,2,:)>=posthoc(2))/n_perm;
    end
    if posthoc(3) < 0
        posthoc_p3 = sum(network_posthoc_null(network_num,1,3,:)<=posthoc(3))/n_perm;
    elseif posthoc(3) > 0
        posthoc_p3 = sum(network_posthoc_null(network_num,1,3,:)>=posthoc(3))/n_perm;
    end
    if posthoc(4) < 0
        posthoc_p4 = sum(network_posthoc_null(network_num,1,4,:)<=posthoc(4))/n_perm;
    elseif posthoc(4) > 0
        posthoc_p4 = sum(network_posthoc_null(network_num,1,4,:)>=posthoc(4))/n_perm;
    end
    network_posthoc_output(k) = strcat(network(network_num),',intra,',num2str(posthoc_p1),',',num2str(posthoc_p2),',',num2str(posthoc_p3),',',num2str(posthoc_p4)); % store permutation p values
    k = k+1;
        
    % [internetwork] calculate observed posthoc (pairwise comparison) values
	pred_val = predict(inter_lme_mdl{network_num},network_lme_tbl);
    posthoc(1) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
    posthoc(2) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
    posthoc(3) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
    posthoc(4) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
		
	% [internetwork] get fraction of permutation values (null distribution) that are as or more extreme than observed posthoc values
    if posthoc(1) < 0
        posthoc_p1 = sum(network_posthoc_null(network_num,2,1,:)<=posthoc(1))/n_perm;
    elseif posthoc(1) > 0
        posthoc_p1 = sum(network_posthoc_null(network_num,2,1,:)>=posthoc(1))/n_perm;
    end
    if posthoc(2) < 0
        posthoc_p2 = sum(network_posthoc_null(network_num,2,2,:)<=posthoc(2))/n_perm;
    elseif posthoc(2) > 0
        posthoc_p2 = sum(network_posthoc_null(network_num,2,2,:)>=posthoc(2))/n_perm;
    end
    if posthoc(3) < 0
        posthoc_p3 = sum(network_posthoc_null(network_num,2,3,:)<=posthoc(3))/n_perm;
    elseif posthoc(3) > 0
        posthoc_p3 = sum(network_posthoc_null(network_num,2,3,:)>=posthoc(3))/n_perm;
    end
    if posthoc(4) < 0
        posthoc_p4 = sum(network_posthoc_null(network_num,2,4,:)<=posthoc(4))/n_perm;
    elseif posthoc(4) > 0
        posthoc_p4 = sum(network_posthoc_null(network_num,2,4,:)>=posthoc(4))/n_perm;
    end
    network_posthoc_output(k) = strcat(network(network_num),',inter,',num2str(posthoc_p1),',',num2str(posthoc_p2),',',num2str(posthoc_p3),',',num2str(posthoc_p4)); % store permutation p values
    k = k+1; 
end

% collect p values for behavioural scores
j = 2;
k = 2;
for behav_num = 1:length(behav)
	
	% get fraction of permutation values (null distribution) that are as or more extreme than observed fixed effect
    idx = find(strcmp(behav_lme_mdl{behav_num}.Coefficients.Name,'timepoint_2')); % for main effect of time
    if behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) < 0
        time_p = sum(behav_time_beta_null(behav_num,:)<=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    elseif behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) > 0
        time_p = sum(behav_time_beta_null(behav_num,:)>=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(behav_lme_mdl{behav_num}.Coefficients.Name,'group_2')); % for main effect of group
    if behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) < 0
        group_p = sum(behav_group_beta_null(behav_num,:)<=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    elseif behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) > 0
        group_p = sum(behav_group_beta_null(behav_num,:)>=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    end
    idx = find(strcmp(behav_lme_mdl{behav_num}.Coefficients.Name,'timepoint_2:group_2')); % for interaction effect
    if behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) < 0
        interaction_p = sum(behav_int_beta_null(behav_num,:)<=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    elseif behav_lme_mdl{behav_num}.Coefficients.Estimate(idx) > 0
        interaction_p = sum(behav_int_beta_null(behav_num,:)>=behav_lme_mdl{behav_num}.Coefficients.Estimate(idx))/n_perm;
    end
    behav_perm_output(j) = strcat(behav(behav_num),',',num2str(time_p),',',num2str(group_p),',',num2str(interaction_p)); % store permutation p values
    j = j+1;
    
    % calculate observed posthoc (pairwise comparison) values
    pred_val = predict(behav_lme_mdl{behav_num},behav_lme_tbl);
    idx_pre = find(behav_lme_tbl.timepoint == '1');
    idx_post = find(behav_lme_tbl.timepoint == '2');
    idx_mind = find(behav_lme_tbl.group == '1'); % MBTI
    idx_slp = find(behav_lme_tbl.group == '2'); % SHEEP
    posthoc(1) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_post,idx_mind))); % pre M - post M
    posthoc(2) = mean(pred_val(intersect(idx_pre,idx_slp)))-mean(pred_val(intersect(idx_post,idx_slp))); % pre S - post S
    posthoc(3) = mean(pred_val(intersect(idx_pre,idx_mind)))-mean(pred_val(intersect(idx_pre,idx_slp))); % pre M - pre S
    posthoc(4) = mean(pred_val(intersect(idx_post,idx_mind)))-mean(pred_val(intersect(idx_post,idx_slp))); % post M - post S
	
	% get fraction of permutation values (null distribution) that are as or more extreme than observed posthoc values
    if posthoc(1) < 0
        posthoc_p1 = sum(behav_posthoc_null(behav_num,1,:)<=posthoc(1))/n_perm;
    elseif posthoc(1) > 0
        posthoc_p1 = sum(behav_posthoc_null(behav_num,1,:)>=posthoc(1))/n_perm;
    end
    if posthoc(2) < 0
        posthoc_p2 = sum(behav_posthoc_null(behav_num,2,:)<=posthoc(2))/n_perm;
    elseif posthoc(2) > 0
        posthoc_p2 = sum(behav_posthoc_null(behav_num,2,:)>=posthoc(2))/n_perm;
    end
    if posthoc(3) < 0
        posthoc_p3 = sum(behav_posthoc_null(behav_num,3,:)<=posthoc(3))/n_perm;
    elseif posthoc(3) > 0
        posthoc_p3 = sum(behav_posthoc_null(behav_num,3,:)>=posthoc(3))/n_perm;
    end
    if posthoc(4) < 0
        posthoc_p4 = sum(behav_posthoc_null(behav_num,4,:)<=posthoc(4))/n_perm;
    elseif posthoc(4) > 0
        posthoc_p4 = sum(behav_posthoc_null(behav_num,4,:)>=posthoc(4))/n_perm;
    end
    behav_posthoc_output(k) = strcat(behav(behav_num),',',num2str(posthoc_p1),',',num2str(posthoc_p2),',',num2str(posthoc_p3),',',num2str(posthoc_p4)); % store permutation p values
    k = k+1;
end

% readable outputs saved as csv files
writetable(cell2table(network_perm_output'),'TO_BE_SPECIFIED.csv');
writetable(cell2table(network_posthoc_output'),'TO_BE_SPECIFIED.csv');
writetable(cell2table(behav_perm_output'),'TO_BE_SPECIFIED.csv');
writetable(cell2table(behav_posthoc_output'),'TO_BE_SPECIFIED.csv');
% null distributions from permutations and permutation indices saved in mat file
save('TO_BE_SPECIFIED.mat','grp_perm_idx','time_perm','behav_time_beta_null','behav_group_beta_null','behav_int_beta_null','intra_time_beta_null','intra_group_beta_null','intra_int_beta_null','inter_time_beta_null','inter_group_beta_null','inter_int_beta_null','behav_posthoc_null','network_posthoc_null');

