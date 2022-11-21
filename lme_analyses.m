% Perform linear mixed modelling (lme) on FC similarity and behavioural measures to get group and timepoint interactions

% This script assumes the FC similarity values, behavioural data and demographic data are ready to use as inputs
% All parts that need to be updated with specific information are indicated with TO_BE_SPECIFIED

clear;clc % clear workspace variables and command window

exclude = {''}; % subjects to be excluded

% get subjects from subject list file
subj_list = 'TO_BE_SPECIFIED.txt';
subj_name = importdata(subj_list);

%% collect behavioural and demographic measures 
subj_data = readtable('TO_BE_SPECIFIED.csv'); % load demographic and behavioural data 
behav = {'TO_BE_SPECIFIED'}; % behavioural scores used

for subj_num = 1:length(subj_name)
    if sum(strcmpi(subj_name(subj_num),exclude)) == 0 % for all subjects other than the one to be excluded
        for timepoint_num = 1:2 % pre and post
            subj_idx(subj_num) = subj_num; % get subject numbers from subject name
            subj_data_idx = find(~cellfun('isempty',strfind(subj_data{:,'Subject'},[subj_name{subj_num} '-' num2str(timepoint_num)]))); % get idx for the subject (one row per subject)
			
			% collect demographic data (age, gender, group) TO_BE_SPECIFIED accordingly based on variable names in demographics file
            demographics(subj_num,1,timepoint_num) = subj_data{subj_data_idx,'age'}; 
            demographics(subj_num,2,timepoint_num) = subj_data{subj_data_idx,'gender'};
            demographics(subj_num,3,timepoint_num) = subj_data{subj_data_idx,'group'};
            
			for behav_num = 1:length(behav) % collect behavioural scores
				all_behav(subj_num,behav_num,timepoint_num) = subj_data{subj_data_idx,behav{behav_num}};
			end			
        end
    else % display message for excluded subject
        exclude_subj = strcat(subj_name(subj_num),': subject excluded');
        disp(exclude_subj);
    end
end

subj_nonzero_idx = find(subj_idx); % get index of subjects included in analysis 

% keep only data for subjects analysed
subj_name = subj_name(subj_nonzero_idx); 
demographics = demographics(subj_nonzero_idx,:,:); 
all_behav = all_behav(subj_nonzero_idx,:,:);

%% lme for behavioural measures

% create tables to use as input for running lme function - add in demographic variables first
behav_lme_tbl = [table([repmat(subj_name,2,1)],'VariableNames',{'subject'}) table([repmat(1,length(subj_name),1);repmat(2,length(subj_name),1)],'VariableNames',{'timepoint'}) array2table([demographics(:,:,1);demographics(:,:,2)],'VariableNames',{'age','gender','group'})];
behav_lme_tbl.group = categorical(behav_lme_tbl.group); 
behav_lme_tbl.gender = categorical(behav_lme_tbl.gender);
behav_lme_tbl.timepoint = categorical(behav_lme_tbl.timepoint);

behav_lme_output = [];
for behav_num = 1:length(behav)
    behav_lme_tbl = [behav_lme_tbl table([all_behav(:,behav_num,1);all_behav(:,behav_num,2)],'VariableNames',behav(behav_num))]; % add behavioural score to table
    behav_lme_mdl{behav_num} = fitlme(behav_lme_tbl,[behav{behav_num} '~age+gender+group*timepoint+(timepoint|subject)']); % run lme 
    behav_lme_output = [behav_lme_output; [table(repmat(behav(behav_num),behav_lme_mdl{behav_num}.NumCoefficients,1)),dataset2table(behav_lme_mdl{behav_num}.Coefficients(:,:))]]; % store outputs from lme
end

%% lme for FC similarity
network = {'TO_BE_SPECIFIED'}; % network names
load('TO_BE_SPECIFIED.mat'); % input file containing FC similarity values

% create tables to use as input for running lme function - add in demographic variables first
network_lme_tbl = [table([repmat(subj_name,2,1)],'VariableNames',{'subject'}) table([repmat(1,length(subj_name),1);repmat(2,length(subj_name),1)],'VariableNames',{'timepoint'}) array2table([demographics(:,:,1);demographics(:,:,2)],'VariableNames',{'age','gender','group'})];
network_lme_tbl.group = categorical(network_lme_tbl.group); 
network_lme_tbl.gender = categorical(network_lme_tbl.gender);
network_lme_tbl.timepoint = categorical(network_lme_tbl.timepoint);

network_lme_output = [];
for network_num = 1:length(network)
	% add FC similarity to table
	network_lme_tbl = [network_lme_tbl table([z_intranetwork_rest_task_sim(:,1,network_num);z_intranetwork_rest_task_sim(:,2,network_num)],'VariableNames',{[network{network_num} '_intra_fc_sim']})]; 
	network_lme_tbl = [network_lme_tbl table([z_internetwork_rest_task_sim(:,1,network_num);z_internetwork_rest_task_sim(:,2,network_num)],'VariableNames',{[network{network_num} '_inter_fc_sim']})];

	% run lme for FC similarity 
	intra_lme_mdl{network_num} = fitlme(network_lme_tbl,[network{network_num} '_intra_fc_sim~age+gender+group*timepoint+(timepoint|subject)']);   
	inter_lme_mdl{network_num} = fitlme(network_lme_tbl,[network{network_num} '_inter_fc_sim~age+gender+group*timepoint+(timepoint|subject)']); 
	
	% store outputs from lme
	network_lme_output = [network_lme_output; table(repmat(network(network_num),intra_lme_mdl{network_num}.NumCoefficients,1),repmat('intra',intra_lme_mdl{network_num}.NumCoefficients,1)),dataset2table(intra_lme_mdl{network_num}.Coefficients(:,:))]; 
	network_lme_output = [network_lme_output; table(repmat(network(network_num),inter_lme_mdl{network_num}.NumCoefficients,1),repmat('inter',inter_lme_mdl{network_num}.NumCoefficients,1)),dataset2table(inter_lme_mdl{network_num}.Coefficients(:,:))]; 
end

%% save outputs

% readable outputs in csv files
writetable(network_lme_output,['TO_BE_SPECIFIED.csv']);
writetable(behav_lme_output,['TO_BE_SPECIFIED.csv']);
% lme models saved in mat file
save(['TO_BE_SPECIFIED.mat'],'behav_lme_mdl','intra_lme_mdl','inter_lme_mdl');
