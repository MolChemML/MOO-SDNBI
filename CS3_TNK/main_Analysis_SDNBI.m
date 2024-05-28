clc; clear;

%--------------------------------------------------------------------------
% Load the simulation
%--------------------------------------------------------------------------
load('CS3_SDNBI_region_EJOR_final.mat');
%load('CS3_SDNBI_210217-dTolsTol_0.005_0.04_mostwell.mat');

PPointsHistory_Re = PPointsHistory;
%--------------------------------------------------------------------------
% Calculate CPU time, HV, CD, DM
%--------------------------------------------------------------------------
for i=1:size(PPointsHistory_Re,1)

    
    [mat_unq,idx_unq]          = unique(round(PPointsHistory_Re(1:i,3:4),5),'rows');
    idx_unq=sortrows(idx_unq,1,'ascend');
    [ndf_index, df_index]      = non_dominated_front(PPointsHistory_Re(idx_unq,3:4)');  
    CS3_Analysis_SDNBI(i,1)   = length(ndf_index);
    CS3_Analysis_SDNBI(i,2:3) = PPointsHistory_Re(i,3:4);          % objective function
    CS3_Analysis_SDNBI(i,4:5) = PPointsHistory_Re(i,5:6);          % decision variables
    CS3_Analysis_SDNBI(i,6:7) = PPointsHistory_Re(i,16:17);        % beta  
    CS3_Analysis_SDNBI(i,8)   = PPointsHistory_Re(i,22);           % dMax 
    CS3_Analysis_SDNBI(i,9:11) = [PPointsHistory_Re(i,20:21),PPointsHistory_Re(i,24)];       % Wall/CPUTime
    CS3_Analysis_SDNBI(i,12)   = sum(PPointsHistory_Re(1:i,20));   % Acc Wall/CPUTime
    CS3_Analysis_SDNBI(i,13)  = sum(PPointsHistory_Re(1:i,21));    % Acc Wall/CPUTime
    
    if i < 3 
       CS3_Analysis_SDNBI(i,14:17) = [0,99,99,99];  %HV, crowding_diststdv,avg,DM
    else
        CS3_Analysis_SDNBI(i,14) = approximate_hypervolume_ms(PPointsHistory_Re(ndf_index,3:4)',[1,1]');
        
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
        
        CS3_Analysis_SDNBI(i,15:17) = [CD_stdv,CD_avg,DM]; 
    end
    CS3_Analysis_SDNBI(i,18)  = sum(PPointsHistory_Re(1:i,24));    % Acc Wall/CPUTime
    CS3_Analysis_SDNBI(i,19)   = PPointsHistory_Re(i,23);           % sMax 
    
end

%save('CS3_mNBI_210217_final_analysis.mat');
save('CS3_SDNBI_region_EJOR_final_analysis.mat');

%hold off
% Plot CPU vs Hypervolume
%plot( CS3_Analysis_SDNBI(3:end,9), CS3_Analysis_SDNBI(3:end,11),'o','MarkerEdgeColor',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',6);
%plot( PF_CS3_SDNBI(:,3), PF_CS3_SDNBI(:,4),'o','MarkerEdgeColor','red','MarkerFaceColor',rgb('Gold'),'LineWidth',1,'MarkerSize',6);
%plot( PF_CS3_SD(:,3), PF_CS3_SD(:,4),'x','MarkerEdgeColor','blue','MarkerFaceColor',rgb('blue'),'LineWidth',1,'MarkerSize',6);
%plot( PF_CS3_mNBI(:,3), PF_CS3_mNBI(:,4),'^','MarkerEdgeColor','green','MarkerFaceColor',rgb('green'),'LineWidth',1,'MarkerSize',6);
