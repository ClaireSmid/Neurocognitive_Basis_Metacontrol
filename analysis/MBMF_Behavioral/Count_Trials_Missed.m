
load('../kidgroupdata.mat');

data = kidgroupdata.subdata;

sums = [];

for i = 1:length(data)
   
    
    pp = data(i);
    
    thisid = pp.id;
    thesepoints = kidgroupdata.points(i);
    lopoints = kidgroupdata.points_bystake(i,1);
    hipoints = kidgroupdata.points_bystake(i,2);
    
    if thisid < 211
        
        rr.id = thisid;
        rr.trials = pp.N;
        rr.missed = sum([pp.timeout]);
        rr.miss_p = sum([pp.timeout])/pp.N;
        rr.Avg_Pts = thesepoints;
        rr.Avg_Pts_lo = lopoints;
        rr.Avg_Pts_hi = hipoints;

        sums = [sums; rr];
        
    end
    
end

T = struct2table(sums);
writetable(T,'MBMF_MissedTrials_31May2022.csv','Delimiter',',')