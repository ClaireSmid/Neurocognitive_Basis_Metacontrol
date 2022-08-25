%% Wrapper for Stay Probabilities
% for both adults, kids and the two models.
% the regression was conducted in R based on the output from these script,
% using the lme4 package. 

clearvars;
clc;
  
% load their combined data, and add some columns 
load('../kidgroupdata.mat');
groupdata_I = make_raw_data(kidgroupdata);

% for the two models
for m = 1:2

    fprintf('Running model no. %d...\n\n', m)

    % for the P6 model
    if m == 1

        % get their stay probabilities
        groupdata = Get_Stay_Probabilities(groupdata_I,m);

        % save their data
        save Stay_Prob_Kids_P6_full.mat groupdata
        writetable(groupdata.table,'../Stay_Prob_Kids_P6.csv','Delimiter',',')
        fprintf('saving data from model %d...\n\n', m)

        % clear it
        clear groupdata

    elseif m == 2

        % get their stay probabilities
        groupdata = Get_Stay_Probabilities(groupdata_I,m);

        % save their data
        save Stay_Prob_Kids_P7_full.mat groupdata
        writetable(groupdata.table,'../Stay_Prob_Kids_P7.csv','Delimiter',',')
        fprintf('saving data from model %d...\n\n', m)

        % clear it
        clear groupdata

    end

end

   