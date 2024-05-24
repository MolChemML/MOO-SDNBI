clc; clear;

%--------------------------------------------------------------------------
% Load the simulation
%--------------------------------------------------------------------------
load('CS5_SDNBI_Run.mat')
PPointsHistory_Re=PPointsHistory;
%--------------------------------------------------------------------------
% Calculate CPU time, HV, CD, DM
%--------------------------------------------------------------------------
for i=1:size(PPointsHistory_Re,1)

    [ndf_index, df_index] = non_dominated_front(PPointsHistory_Re(1:i,3:4)');
    
    CS5_Analysis_SDNBI(i,1)   = length(ndf_index);
    CS5_Analysis_SDNBI(i,2:3) = PPointsHistory_Re(i,3:4);          % objective function
    CS5_Analysis_SDNBI(i,4:5) = PPointsHistory_Re(i,5:6);          % decision variables
    CS5_Analysis_SDNBI(i,6:7) = PPointsHistory_Re(i,21:22);        % beta  
    CS5_Analysis_SDNBI(i,8)   = PPointsHistory_Re(i,end-2);           % dMax 
    CS5_Analysis_SDNBI(i,9:11) = [PPointsHistory_Re(i,end-4:end-3),PPointsHistory_Re(i,end)];       % Wall/CPUTime
    CS5_Analysis_SDNBI(i,12)   = sum(PPointsHistory_Re(1:i,end-4));   % Acc Wall/CPUTime
    CS5_Analysis_SDNBI(i,13)  = sum(PPointsHistory_Re(1:i,end-3));    % Acc Wall/CPUTime
    
    if i < 3 
       CS5_Analysis_SDNBI(i,14:17) = [0,99,99,99];  %HV, crowding_diststdv,avg,DM
    else
        CS5_Analysis_SDNBI(i,14) = approximate_hypervolume_ms(PPointsHistory_Re(ndf_index,3:4)',[1,1]');
        
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
        
        CS5_Analysis_SDNBI(i,15:17) = [CD_stdv,CD_avg,DM]; 
    end
    CS5_Analysis_SDNBI(i,18)  = sum(PPointsHistory_Re(1:i,end));    % Acc Wall/CPUTime
    CS5_Analysis_SDNBI(i,19)  = PPointsHistory_Re(i,end-1);    % Acc Wall/CPUTime 
end


save('CS5_SDNBI_final_analysis.mat');

hold off

plot( CS5_Analysis_SDNBI(3:end,9), CS5_Analysis_SDNBI(3:end,11),'o','MarkerEdgeColor',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',6);
