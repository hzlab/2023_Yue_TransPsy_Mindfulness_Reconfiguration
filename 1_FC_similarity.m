% Calculate FC similarity between rest and task FC

% This script assumes the FC matrices for rest and task are already generated and ready to use as inputs
% With individual files for rest and task FC matrices from each subject at each timepoint
% All parts that need to be updated with specific information are indicated with TO_BE_SPECIFIED

clear;clc % clear workspace variables and command window

exclude = {''}; % subjects to be excluded

output_folder = 'TO_BE_SPECIFIED';
mkdir(output_folder)
output_file = 'TO_BE_SPECIFIED.mat';

% folders for subject resting and task FC data
rest_folder = {'TO_BE_SPECIFIED'};
task_folder = {'TO_BE_SPECIFIED'};

% get subjects from subject list file
subj_list = 'TO_BE_SPECIFIED.txt';
subj_name = importdata(subj_list);

network = {'TO_BE_SPECIFIED'}; % network names
network_roi = load('TO_BE_SPECIFIED.mat'); % ROI assignment to each network
max_roi = TO_BE_SPECIFIED'; % this should be the total number of ROIs 

% Preallocate some variables 
intranetwork_rest_task_sim = NaN(length(subj_name),2,length(network)); % (subject,timepoint,network)
internetwork_rest_task_sim = NaN(length(subj_name),2,length(network)); % (subject,timepoint,network)
subj_idx = zeros(length(subj_name),1); % (subject)

% generate rest and task network-based FC matrices from ROI-based FC matrices, then calculate rest-task FC similarity
for subj_num = 1:length(subj_name)
    if sum(strcmpi(subj_name(subj_num),exclude)) == 0 % for all subjects other than the one to be excluded
	
        for timepoint_num = 1:2 % pre and post
            % filenames for rest and task fc matrices (ROI-based) for each subject
            subj_rest_fc = strcat(rest_folder,subj_name(subj_num),'-',num2str(timepoint_num),'/restfc_mats.mat'); % TO_BE_SPECIFIED accordingly based on folder structure and file name
            subj_task_fc = strcat(task_folder,subj_name(subj_num),'-',num2str(timepoint_num),'/taskfc_mats.mat'); % TO_BE_SPECIFIED accordingly based on folder structure and file name

            if exist(subj_rest_fc{:}, 'file') && exist(subj_task_fc{:}, 'file') % check that rest and task FC matrices exist
                subj_idx(subj_num) = subj_num; 
                
                % get rest and task FC matrices (ROI-based, fisher's z transformed) for each subject
                load(subj_rest_fc{:});
                rest_fc = z_fc_mat; % TO_BE_SPECIFIED accordingly based on variable name
                load(subj_task_fc{:});
                task_fc = z_fc_mat; % TO_BE_SPECIFIED accordingly based on variable name
                
                for network_num = 1:length(network)

                    % calculate intranetwork FC (only need values above main diagonal --> symmetrical matrix + exclude self-connections of ROI with itself)
                    mask = triu(ones(length(network_roi{network_num})),1); % generate mask that has 1 above main diagonal and 0 elsewhere
                    mask(mask==0) = NaN; % set 0s to NaNs (everywhere except above main diagonal)
                    
					% apply mask to intranetwork ROI-pairs from FC matrices --> get only values above main diagonal from FC matrices and change the rest to NaN
                    intranetwork_rest_fc = reshape(mask.*rest_fc(network_roi{network_num},network_roi{network_num}),length(network_roi{network_num})*length(network_roi{network_num}),1); % vectorized
                    intranetwork_task_fc = reshape(mask.*task_fc(network_roi{network_num},network_roi{network_num}),length(network_roi{network_num})*length(network_roi{network_num}),1); % vectorized
					intranetwork_rest_task_sim(subj_num,timepoint_num,network_num) = corr(intranetwork_rest_fc,intranetwork_task_fc,'rows','complete');
										
                    % extract the relevant internetwork ROI-pairs from FC matrices --> get whole row in FC matrices
                    inter_idx = setdiff(1:max_roi,network_roi{network_num}); 
                    internetwork_rest_fc = reshape(rest_fc(network_roi{network_num},inter_idx),length(network_roi{network_num})*length(inter_idx),1); % vectorized
                    internetwork_task_fc = reshape(task_fc(network_roi{network_num},inter_idx),length(network_roi{network_num})*length(inter_idx),1); % vectorized
                    clear inter_idx
                    internetwork_rest_task_sim(subj_num,timepoint_num,network_num) = corr(internetwork_rest_fc,internetwork_task_fc,'rows','complete');
                    
                    clear intranetwork_rest_fc intranetwork_task_fc internetwork_rest_fc internetwork_task_fc
                end
                
                clear rest_fc task_fc
                
            else % rest and/or task matrices do not exist for the subject --> display message
                missing_subj_fc = strcat(subj_name(subj_num), ': missing FC matrix/matrices');
                disp(missing_subj_fc);
            end
        end
    else % display message for excluded subject
        exclude_subj = strcat(subj_name(subj_num),': subject excluded');
        disp(exclude_subj);
    end
end

subj_nonzero_idx = find(subj_idx); % get index of subjects included in analysis (out of all subjects with good motion qc)
subj_name = subj_name(subj_nonzero_idx); % keep only names/identifiers for subjects analysed
subj_idx = nonzeros(subj_idx); % remove zeros (for subjects with missing FC/excluded)
intranetwork_rest_task_sim = reshape(intranetwork_rest_task_sim(~isnan(intranetwork_rest_task_sim)),length(subj_idx),2,length(network)); % remove NaN values (subjects with missing FC/excluded)
internetwork_rest_task_sim = reshape(internetwork_rest_task_sim(~isnan(internetwork_rest_task_sim)),length(subj_idx),2,length(network)); % remove NaN values (subjects with missing FC/excluded)

% Fisher's z transformation of rest-task FC correlation
z_intranetwork_rest_task_sim = 0.5.*log((1+intranetwork_rest_task_sim)./(1-intranetwork_rest_task_sim)); %(subject,load,network)
z_internetwork_rest_task_sim = 0.5.*log((1+internetwork_rest_task_sim)./(1-internetwork_rest_task_sim)); %(subject,load,network)

save(output_file,'z_intranetwork_rest_task_sim','z_internetwork_rest_task_sim','subj_name');
