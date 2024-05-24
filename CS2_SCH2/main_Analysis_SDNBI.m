clc; clear;

%--------------------------------------------------------------------------
% Load the simulation
%--------------------------------------------------------------------------
load('CS2_SDNBI_final_Run.mat');
PPointsHistory_Re = PPointsHistory;
CS2_Analysis_SDNBI=[];
%--------------------------------------------------------------------------
% Calculate CPU time, HV, CD, DM
%--------------------------------------------------------------------------
for i=1:size(PPointsHistory_Re,1)
 
    [mat_unq,idx_unq]          = unique(round(PPointsHistory_Re(1:i,3:4),5),'rows');
    [ndf_index, df_index]      = non_dominated_front(PPointsHistory_Re(idx_unq,3:4)');
    CS2_Analysis_SDNBI(i,1)    = length(ndf_index);
    CS2_Analysis_SDNBI(i,2:3)  = PPointsHistory_Re(i,3:4);          % objective function
    CS2_Analysis_SDNBI(i,4:5)  = PPointsHistory_Re(i,5:6);          % decision variables
    CS2_Analysis_SDNBI(i,6:7)  = PPointsHistory_Re(i,16:17);        % beta  
    CS2_Analysis_SDNBI(i,8)    = PPointsHistory_Re(i,22);           % dMax 
    CS2_Analysis_SDNBI(i,9:11) = [PPointsHistory_Re(i,20:21),PPointsHistory_Re(i,24)];       % Wall/CPUTime
    CS2_Analysis_SDNBI(i,12)   = sum(PPointsHistory_Re(1:i,20));   % Acc Wall/CPUTime
    CS2_Analysis_SDNBI(i,13)   = sum(PPointsHistory_Re(1:i,21));    % Acc Wall/CPUTime
    
    if i < 3 
       CS2_Analysis_SDNBI(i,14:17) = [0,99,99,99];  %HV, crowding_diststdv,avg,DM
    else
        CS2_Analysis_SDNBI(i,14) = approximate_hypervolume_ms(PPointsHistory_Re(ndf_index,3:4)',[1,1]');
        
        crowding_dist=compute_crowding_distances(PPointsHistory_Re(ndf_index,3:4)');
        CD = [];
        for j=1:length(crowding_dist)
            if crowding_dist(j) > 1e5;
                continue
            end
            CD(1,end+1) = crowding_dist(j);
        end
        CD_stdv=std(CD');
        CD_avg =mean(CD');
        DM = compute_quality_pareto(PPointsHistory_Re(ndf_index,3:4), [1,1], [0,0]);
        
        CS2_Analysis_SDNBI(i,15:17) = [CD_stdv,CD_avg,DM]; 
    end
    CS2_Analysis_SDNBI(i,18)  = sum(PPointsHistory_Re(1:i,24));    % Acc Wall/CPUTime
    
end


save('CS2_SDNBI_final_analysis.mat');

n_fixed=100;
[unq,id_unq,dummy] = unique(round(PPointsHistory_Re(1:n_fixed,3:4),6),'rows');
 id_unq = sortrows(id_unq,1,'ascend');
[ndf,df]= non_dominated_front(PPointsHistory_Re(id_unq,3:4)');
 PF_fixIter_CS2_SDNBI = PPointsHistory_Re(ndf,3:4);

