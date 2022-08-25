%% Notes for function
% Claire Smid, August 2020
% This function starts a parfor loop to fit participants simultaneously. 

function [kidresults] = Fit_model_wrapper_parfor_P6

    
load kidgroupdata.mat
participant_data = kidgroupdata;

%% Step 1. Start analysing all subjects with a parfor loop
subs = length(participant_data.subdata(:));

% random starting locations for optimizer
nstarts = 100; % 100;

parfor i = 1:subs

    data = participant_data.subdata(i);
    data.rews = data.rews.*9;
    fprintf('Fitting participant no. %d of sample %d\n\n',i)

    % run optimization
    params = set_params_P6;
    f = @(x,data) MB_MF_llik_P6(x,data);
    results(i) = mfit_optimize(f,params,data,nstarts);

end

%% Step 2. Save results per sample
kidresults = results;
save kidsresults_P6_T0_5April2022.mat kidresults


end
   