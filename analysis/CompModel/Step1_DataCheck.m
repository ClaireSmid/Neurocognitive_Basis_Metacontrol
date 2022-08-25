%%% first script to create the combined data files.

% 1. combine data files into one matrix
% exclude any participant with less than 90 trials and >30% missing trials(2sd from the mean missing trials)

addpath('data')


disp('Step 1. Combining data files.\n\n\n')
% need to create a log of removed participants here.
[kidgroupdata, kidlog] = groupanalysis;


function [kidgroupdata, kidlog] = groupanalysis

    
groupdata = [];
com_log = [];

files = dir('data/ModelBasedLearning_s*');
disp(files)
kidgroupdata = [];

nrfiles = length(files);

s = 0;

for i = 1:nrfiles
   
    fprintf('Analysing PP %d\n',i)
    
    load([files(i).folder,'/',files(i).name]);
    
    subdata = [];
    log = [];
    log.id = str2double(data.id);
    log.trials_attempted = sum(data.block>0);
    
    if sum(data.block>0) >= 90
        
        log.init_inclusion = 1;
        
        trials = 1:sum(data.block>0);
        
        % general info
        subdata.id = str2double(data.id);
        subdata.age = data.age;
        subdata.data = data.date;
        subdata.dataFileName = data.dataFileName;
        
        % task info
        subdata.rews = data.rews(trials,:)/9;
        subdata.block = data.block(trials);
        subdata.s = data.s(trials,:);
        subdata.stimuli = data.stimuli(trials,:);
        subdata.stake = data.stake(trials);
        subdata.choice = data.choice(trials);
        subdata.rt = data.rt(trials,:);
        subdata.score = data.score(trials)/9;
        subdata.points = data.points(trials)/9;
        subdata.timeout = data.timeout(trials,:);
        subdata.rocketOrder = data.rocketOrder;
        subdata.purpleOrder = data.purpleOrder;
        subdata.redOrder = data.redOrder;
        
        subdata.N = length(subdata.block);
        subdata.missed = any(subdata.timeout,2);
        
        if sum(subdata.timeout(:)) > 0
            Missers = sum(subdata.timeout(:));
            P_missed = (Missers*100)/max(trials);
        else
            P_missed = 0;
        end
        
        log.P_missed = P_missed;
        log.actual_trials = (sum(data.block>0) - sum(subdata.missed));

        if P_missed < 30 % so if they did 140 trials and missed less than 42 trials, they will be included
%         if mean(subdata.missed) < 0.3 % throw out participants with > 30% missed data
            
            log.mean_missed_trials = mean(subdata.missed);
            log.missed = sum(subdata.timeout(:,1));
%             log.scnd_stage_missed = sum(subdata.timeout(:,2));
            log.final_inclusion = 1;
            
            % included subject
            s = s + 1;
            
            groupdata.subdata(s) = subdata;
%             id_num(s) = str2double(data.id);
            
        else
            
            log.mean_missed_trials = mean(subdata.missed);
            log.missed = sum(subdata.timeout(:,1));
%             log.scnd_stage_missed = sum(subdata.timeout(:,2));
            log.final_inclusion = 0;
            
        end
        
    else
        log.init_inclusion = 0;
        log.P_missed = NaN;
        log.actual_trials = NaN;
        log.mean_missed_trials = NaN;
        log.missed = NaN;
%         log.scnd_stage_missed = NaN;
        log.final_inclusion = 0;
        
    end
    
    % save log
    com_log = [com_log, log];
    
    
    subdata = groupdata.subdata(s);
    
    % save some info per pp
    groupdata.N(s,1) = subdata.N;
    groupdata.missed(s,1) = mean(subdata.missed);
    groupdata.age(s,1) = subdata.age;
    groupdata.id(s,1) = subdata.id;
    
    % calculate corrected reward rate here
    groupdata.points(s,1) = mean(subdata.points(~subdata.missed)) - mean(mean(subdata.rews(~subdata.missed,:)));
    groupdata.points_bystake(s,1) = mean(subdata.points(subdata.stake==1&~subdata.missed)) - mean(mean(subdata.rews(subdata.stake==1&~subdata.missed,:)));
    groupdata.points_bystake(s,2) = mean(subdata.points(subdata.stake==5&~subdata.missed)) - mean(mean(subdata.rews(subdata.stake==5&~subdata.missed,:)));
    
    % below lines are just the same as above but done separately to check. 
%     groupdata.rewardrate_rr(s,1) = mean(subdata.points(~subdata.missed)); % these are the points not counting trials where they made no decision
%     groupdata.avg_rew(s,1) = mean(mean(subdata.rews(~subdata.missed,:)));
%     groupdata.corr_rr(s,1) = groupdata.rewardrate_rr(s,1) - groupdata.avg_rew(s,1);
%     groupdata.rr_bystake(s,1) = mean(subdata.points(subdata.stake==1&~subdata.missed));
%     groupdata.rr_bystake(s,2) = mean(subdata.points(subdata.stake==5&~subdata.missed));
%     groupdata.avg_rew_bystake(s,1) = mean(mean(subdata.rews(subdata.stake==1&~subdata.missed,:)));
%     groupdata.avg_rew_bystake(s,2) = mean(mean(subdata.rews(subdata.stake==5&~subdata.missed,:)));
%     groupdata.corr_rr_bystake(s,1) = groupdata.rr_bystake(s,1) - groupdata.avg_rew_bystake(s,1);
%     groupdata.corr_rr_bystake(s,2) = groupdata.rr_bystake(s,2) - groupdata.avg_rew_bystake(s,2);
    
end

% save and rename here
% all_included = sum([com_log.final_inclusion]);
    
% disp('How many pps included?: ');
% disp(all_included);


kidgroupdata = groupdata;
save kidgroupdata.mat kidgroupdata
kidlog = com_log;
save kidlog.mat kidlog


end