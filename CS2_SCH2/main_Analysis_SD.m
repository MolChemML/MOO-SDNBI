clc; clear;

%--------------------------------------------------------------------------
% Load the simulation
%--------------------------------------------------------------------------
load('CS2_SD_final_run.mat');

%--------------------------------------------------------------------------
% Calculate CPU time, HV, CD, DM
%--------------------------------------------------------------------------
for i=1:size(PPoints_History,1)
    [mat_unq,idx_unq]=unique(round(PPoints_History(1:i,3:4),5),'rows');
    [ndf_index, df_index] = non_dominated_front(PPoints_History(idx_unq,3:4)');
    CS2_Analysis_SD(i,1)   = length(ndf_index);
    CS2_Analysis_SD(i,2:3) = PPoints_History(i,3:4);          % objective function
    CS2_Analysis_SD(i,4:5) = PPoints_History(i,5:6);          % decision variables
    CS2_Analysis_SD(i,6:7) = PPoints_History(i,1:2);          % weight  
    CS2_Analysis_SD(i,8)   = PPoints_History(i,end-1);         % dMax 
    CS2_Analysis_SD(i,9:11) = [PPoints_History(i,12:13),PPoints_History(i,end)];    % Wall/CPUTime
    CS2_Analysis_SD(i,12)   = sum(PPoints_History(1:i,12));   % Acc Wall/CPUTime
    CS2_Analysis_SD(i,13)  = sum(PPoints_History(1:i,13));    % Acc Wall/CPUTime
    
    if i < 3 
       CS2_Analysis_SD(i,14:17) = [0,99,99,99];  %HV, crowding_diststdv,avg,DM
    else
        CS2_Analysis_SD(i,14) = approximate_hypervolume_ms(PPoints_History(ndf_index,3:4)',[1,1]');
        
        crowding_dist=compute_crowding_distances(PPoints_History(ndf_index,3:4)');
        CD = [];
        for j=1:length(crowding_dist)
            if crowding_dist(j) > 1e5;
                continue
            end
            CD(1,end+1) = crowding_dist(j);
        end
        CD_stdv=std(CD');
        CD_avg =mean(CD');
        DM = compute_quality_pareto(PPoints_History(ndf_index,3:4), [1,1], [0,0]);
        
        CS2_Analysis_SD(i,15:17) = [CD_stdv,CD_avg,DM]; 
    end

    CS2_Analysis_SD(i,18)  = sum(PPoints_History(1:i,15));    % Acc Wall/CPUTime
    
end

save('CS2_SD_final_analysis.mat');

