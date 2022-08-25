function groupdata = make_raw_data(data)

groupdata = [];


for i = 1:length(data.subdata)

    thisdata = data.subdata(i);

    groupdata.subdata(i).id = thisdata.id(1);
    groupdata.subdata(i).state1 = thisdata.s(:,1);
    groupdata.subdata(i).stim_left = thisdata.stimuli(:,1);
    groupdata.subdata(i).stim_right = thisdata.stimuli(:,2);
    groupdata.subdata(i).rt = thisdata.rt(:,1);
    groupdata.subdata(i).choice = thisdata.choice(:);
%     groupdata.subdata(s).response1 = subdata(:,response1_i);
%     groupdata.subdata(i).rt2 = thisdata.rt(:,2);
    groupdata.subdata(i).rew1 = thisdata.rews(:,1);
    groupdata.subdata(i).rew2 = thisdata.rews(:,2);
    groupdata.subdata(i).rewdiff = abs(thisdata.rews(:,1) - thisdata.rews(:,2));
    groupdata.subdata(i).points = thisdata.points(:);
    groupdata.subdata(i).state2 = thisdata.s(:,2);
    groupdata.subdata(i).stake = thisdata.stake(:);
    groupdata.subdata(i).score = thisdata.score(:);
%     groupdata.subdata(s).practice = subdata(:,practice_i);
    groupdata.subdata(i).rews(:,1) = thisdata.rews(:,1);
    groupdata.subdata(i).rews(:,2) = thisdata.rews(:,2);
    groupdata.subdata(i).trial = [1:length(thisdata.choice)]';
    
    % calculate for previous trials, for the logistic regression
    % (e.g. points achieved on t-1 trial to predict staying on trial t)
    groupdata.subdata(i).missed = (thisdata.rt(:,1) == -1);
    groupdata.subdata(i).prevmissed = [1; groupdata.subdata(i).missed(1:(length(thisdata.choice)-1))];
    groupdata.subdata(i).prevchoice1 = [0; groupdata.subdata(i).choice(1:(length(groupdata.subdata(i).choice)-1))];
    groupdata.subdata(i).prevstate1 = [0; groupdata.subdata(i).state1(1:(length(groupdata.subdata(i).state1)-1))];
    groupdata.subdata(i).prevstate2 = [0; groupdata.subdata(i).state2(1:(length(groupdata.subdata(i).state2)-1))];
    groupdata.subdata(i).prevstake = [0; groupdata.subdata(i).stake(1:(length(groupdata.subdata(i).stake)-1))];
    groupdata.subdata(i).prevpoints = [0; groupdata.subdata(i).points(1:(length(groupdata.subdata(i).points)-1))];
    groupdata.subdata(i).prevrewdiff = [0; groupdata.subdata(i).rewdiff(1:(length(groupdata.subdata(i).rewdiff)-1))];
    
    groupdata.subdata(i).same = double(groupdata.subdata(i).prevstate1 == groupdata.subdata(i).state1);
    groupdata.subdata(i).stay = double(groupdata.subdata(i).prevstate2 == groupdata.subdata(i).state2); % based on planet choice
    
end

end
