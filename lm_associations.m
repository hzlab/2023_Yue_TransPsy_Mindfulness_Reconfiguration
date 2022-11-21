% use linear models (lm) to look at relationship between changes in behaviour and changes in FC similarity

% This script assumes the lme models have been run and saved for FC similarity and behavioural scores to use as inputs

clear;clc % clear workspace variables and command window

with_group = 1; % include group term (1) in lm or not (0)

load('TO_BE_SPECIFIED.mat') % input file containing FC similarity values
load('TO_BE_SPECIFIED.mat') % input file containing lme models
clear all_behav

network = {'TO_BE_SPECIFIED'}; % network names
behav = {'TO_BE_SPECIFIED'}; % behavioural scores used

% get behavioural change
tmp_behav = behav_lme_mdl{length(behav)}.Variables;
all_behav(:,:,1) = table2array(tmp_behav(tmp_behav.timepoint=='1',6:end)); % only take behavioural scores from table
all_behav(:,:,2) = table2array(tmp_behav(tmp_behav.timepoint=='2',6:end)); % only take behavioural scores from table
behav_lm_tbl = [];
for behav_num = 1:length(behav)
   behav_lm_tbl = [behav_lm_tbl table(all_behav(:,behav_num,2)-all_behav(:,behav_num,1),'VariableNames',strcat(behav(behav_num),'_diff'))]; % difference in behavioural score
end

% add demographics and behavioural change to table for lm
network_lm_tbl = [table(subj_name,'VariableNames',{'subject'}) array2table(demographics(:,:,1),'VariableNames',{'age','gender','group'}) behav_lm_tbl]; 
network_lm_tbl.group = categorical(network_lm_tbl.group); 
network_lm_tbl.gender = categorical(network_lm_tbl.gender);

% get FC similarity change and run lm 
network_lm_output = [];
for network_num = 1:length(network)
	% add FC similarity change to the table
	network_lm_tbl = [network_lm_tbl table([z_intranetwork_rest_task_sim(:,2,network_num)-z_intranetwork_rest_task_sim(:,1,network_num)],'VariableNames',{[network{network_num} '_intra_fc_sim_diff']}) table([z_internetwork_rest_task_sim(:,2,network_num)-z_internetwork_rest_task_sim(:,1,network_num)],'VariableNames',{[network{network_num} '_inter_fc_sim_diff']})];
    for behav_num = 1:length(behav)
		% run lm for associations
		if with_group == 1
			intra_lm_mdl{network_num,behav_num} = fitlm(network_lm_tbl,[behav{behav_num} '_diff~age+gender+group*' network{network_num} '_intra_fc_sim_diff']);
			inter_lm_mdl{network_num,behav_num} = fitlm(network_lm_tbl,[behav{behav_num} '_diff~age+gender+group*' network{network_num} '_inter_fc_sim_diff']);
		else
			intra_lm_mdl{network_num,behav_num} = fitlm(network_lm_tbl,[behav{behav_num} '_diff~age+gender+' network{network_num} '_intra_fc_sim_diff']);
			inter_lm_mdl{network_num,behav_num} = fitlm(network_lm_tbl,[behav{behav_num} '_diff~age+gender+' network{network_num} '_inter_fc_sim_diff']);
		end
		% collect output from lm
        network_lm_output = [network_lm_output; [cell2table(repmat({'intra'},size(intra_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var1'}),cell2table(repmat(network(network_num),size(intra_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var2'}),cell2table(repmat(behav(behav_num),size(intra_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var3'}),intra_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,intra_lm_mdl{network_num,behav_num}.Coefficients]];
        network_lm_output.Properties.RowNames = {};
        network_lm_output = [network_lm_output; [cell2table(repmat({'inter'},size(inter_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var1'}),cell2table(repmat(network(network_num),size(inter_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var2'}),cell2table(repmat(behav(behav_num),size(inter_lm_mdl{network_num,behav_num}.Coefficients,1),1),'VariableNames',{'Var3'}),inter_lm_mdl{network_num,behav_num}.Coefficients.Properties.RowNames,inter_lm_mdl{network_num,behav_num}.Coefficients]];
        network_lm_output.Properties.RowNames = {};
    end
end

if with_group == 1
	writetable(network_lm_output,'TO_BE_SPECIFIED.csv'); % save readable output from lm
	save('TO_BE_SPECIFIED.mat','intra_lm_mdl','inter_lm_mdl'); % save linear models 
else
	writetable(network_lm_output,'TO_BE_SPECIFIED.csv'); % save readable output from lm
	save('TO_BE_SPECIFIED.mat','intra_lm_mdl','inter_lm_mdl'); % save linear models 
end
