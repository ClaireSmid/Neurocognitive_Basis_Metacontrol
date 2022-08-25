%% Stay Probability script

function groupdata = Get_Stay_Probabilities(data,m)
% this is based on Wouter's script from 2016

nrsubs = length(data.subdata);

groupdata.subdata = data.subdata;

T = [];
C = [];


for s = 1:nrsubs
    
    thisdata = groupdata.subdata(s);
    groupdata.id(s) = thisdata.id;
    
    fprintf('Analysing PP no. %d, pp code: %d\n',s,groupdata.id(s))

    % remove missing trials
    stake = thisdata.stake(~(thisdata.missed|thisdata.prevmissed));
    same = thisdata.same(~(thisdata.missed|thisdata.prevmissed));
    stay = thisdata.stay(~(thisdata.missed|thisdata.prevmissed));
    prevpoints = thisdata.prevpoints(~(thisdata.missed|thisdata.prevmissed));
    prevrewdiff = thisdata.prevrewdiff(~(thisdata.missed|thisdata.prevmissed));

   
    stake(stake==1) = -1;
    stake(stake==5) = 1; % added
    same = double(same);
    same(same==0) = -1;
    prevpoints = double(prevpoints);
    prevrewdiff = double(prevrewdiff);

    % either include stake or not, for the 6 and 7-parameter models
    if m == 1
        T = [T; groupdata.id(s)*ones(length(same),1) prevrewdiff prevpoints same stay];
    elseif m == 2
        T = [T; groupdata.id(s)*ones(length(same),1) prevrewdiff prevpoints same stake stay];
    end
    
   
end

% combine the data
if m == 1 % no stake 
    groupdata.table = table(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),'VariableNames',{'subnr' 'prevrewdiff' 'prevpoints' 'same' 'stay'});

elseif m == 2 % yes stake
    groupdata.table = table(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),T(:,6),'VariableNames',{'subnr' 'prevrewdiff' 'prevpoints' 'same' 'stake' 'stay'});
end



end

