clear; clc;

%--------------------------------------------------------------------------
% Load the simulation
%--------------------------------------------------------------------------

load('CS2_mNBI_final_run.mat');

%--------------------------------------------------------------------------
% Calculate CPU time, HV, CD, DM
%--------------------------------------------------------------------------
for i=1:size(PPoints_History,1)
 
    [mat_unq,idx_unq]          = unique(round(PPoints_History(1:i,3:4),5),'rows');
    [ndf_index, df_index]      = non_dominated_front(PPoints_History(idx_unq,3:4)');  
    CS2_Analysis_mNBI(i,1)    = length(ndf_index);
    CS2_Analysis_mNBI(i,2:3) = PPoints_History(i,3:4);          % objective function
    CS2_Analysis_mNBI(i,4:5) = PPoints_History(i,5:6);          % decision variables
    CS2_Analysis_mNBI(i,6:7) = PPoints_History(i,16:17);        % beta  
    CS2_Analysis_mNBI(i,8:9) = PPoints_History(i,end-1:end);    % Wall/CPUTime
    CS2_Analysis_mNBI(i,10)   = sum(PPoints_History(1:i,end-1)); % Wall/CPUTime
    CS2_Analysis_mNBI(i,11)  = sum(PPoints_History(1:i,end));   % Wall/CPUTime
    
    if i < 3 
       CS2_Analysis_mNBI(i,12:15) = [0,99,99,99];  %HV, crowding_diststdv,avg,DM
    else
        CS2_Analysis_mNBI(i,12) = approximate_hypervolume_ms(PPoints_History(ndf_index,3:4)',[1,1]');
        
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
        
        CS2_Analysis_mNBI(i,13:15) = [CD_stdv,CD_avg,DM]; 
    end

    
end

save('CS2_mNBI_final_analysis.mat');


        