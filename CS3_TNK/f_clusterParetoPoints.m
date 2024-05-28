function [data_cluster,b_test, data_outer] = f_clusterParetoPoints(index_sub,PPoints,A,b,z_star,dim, c_outer)
    % Cluster the pareto points
    % Only consider the points that define the current facet, no combining
    % proecdure with adjacent facets, not include t in outer approximation
    
    data_cluster = {}; data_outer = {}; i_devide2=1;
    IsConvex = index_sub(1);  % -3: convex , -9: noncovex (concave) 
    multiplier = (IsConvex == -9) * -1 + (IsConvex ~= -9); % Set multiplier based on convexity
    
    flag_extreme = getFlagExtreme(index_sub);
    
    % Extract points and normals
    z = PPoints(index_sub(2:end),dim+1:2*dim); % Current Pareto points by index
    z = [z;z_star];
    a_normal = [A(index_sub(2:end),:); A(end,:)]; % Select corresponding facet normals
    b_offset = [b(index_sub(2:end),:); b(end,:)]; % Select corresponding facet offsets

    
    data_cluster = {[0 0]}; data_outer = {[0 0]};
    [z_sorted,k_boundary]=sortrows(z,1,'ascend');
   
    % Iterate through points to form clusters (always at least triplet)
    for n=2:size(z_sorted,1)-1
       tPoints =  [z_sorted(n-1,:);z_sorted(n,:);z_sorted(n+1,:)];      
       flag_subdetect = 0; b_test2=[];
       
       % Check convexity/concavity for each subset
       for j=1:size(tPoints,1)
           iPoints =  tPoints;
           iPoints(j,:) = [];

           pt_test = sortrows(iPoints,1,'ascend'); %pt_test = sortrows(z,1,'ascend');
       
           w_correct = abs(a_normal(k_boundary(n+j-2),:));   
           b2  = abs(b_offset(k_boundary(n+j-2),:));
           w_test = repmat(w_correct,size(pt_test,1),1);
      % b_test = sum(w_test.*pt_test,2) - repmat(multiplier*b,size(pt_test,1),1);
           b_test2(end+1,:) = transpose(sum(w_test.*pt_test,2) - repmat(b2,size(pt_test,1),1));

           if j==2
               b_test = b_test2(j,:)';
           end
    
           c_devide = prod(b_test2(end,:));
           if c_devide < 0
               b_test = b_test2(j,:)';
               flag_subdetect=1;
               break;
           end  
       end
       
       % Handle extreme points adjustment
       b_test2 = adjustExtremePoints(b_test2, flag_extreme, n, size(z_sorted, 1));
       
       % Handle sub-clusters 
       test_onemore = (b_test2(:,1)' > 0 );
       if ~(nnz(test_onemore)== length(test_onemore) || nnz(test_onemore)== 0)
           c_devide = -1; 
           flag_subdetect=1;
           [row,col]  = find(test_onemore > 0);
           b_test = b_test2(col(1),:)';
       end
       
       % Append or devide clusters
       k_boundary2 = []; 
       if c_devide > 0 && flag_subdetect==0;     % all are interior points of the outer approximation 
          % Assign points in the cluster
          if b_test(1,end) > 0
             index_convex = -3;   % convex index
          else
             index_convex=-9;     % nonconex index
          end
        
        for k=i_devide2:size(k_boundary,1)  % matching z to original index in PPoints 
            i_match = find( round(z(k_boundary(k,1),1),6) == round(PPoints(:,dim+1),6) );
            if length(i_match) > 1
              i_match = find( round(z(k_boundary(k,1),2),6) == round(PPoints(:,dim+2),6) );
            end
            k_boundary2(end+1,1)=i_match;
        end
          data_cluster{end,1}=[index_convex;k_boundary2]'; % data_cluster{end+1,1}=[index_convex;k_boundary2]'; 
          data_outer{end,1}=[index_convex;k_boundary2]'; % data_outer{end+1,1}=[index_convex;k_boundary2]'; 

      elseif c_devide < 0 && flag_subdetect==1;       % cluster should be devided into two new clusters since the sign is opposite
          % Devide and recluster the points 
          if b_test(1,end) > 0
              index_convex1 = -3;
              index_convex2 = -9;
           else
             index_convex1 = -9;
             index_convex2 = -3;
           end
           for k=i_devide2:size(k_boundary,1)  % matching z to original index in PPoints 
              i_match = find( round(z(k_boundary(k,1),1),6) == round(PPoints(:,dim+1),6) );
              if length(i_match) > 1
              i_match = find( round(z(k_boundary(k,1),2),6) == round(PPoints(:,dim+2),6) );
              end
              k_boundary2(end+1,1)=i_match;
           end       
            i_devide2 = n ; i_devide = find(round(z_sorted(n,1),6)==round(PPoints(k_boundary2,dim+1),6)); %i_devide = 2 ;
           sub_01=k_boundary2(1:i_devide,:); sub_02=k_boundary2(i_devide:end,:);
           data_cluster{end,1}=[index_convex1; sub_01]';% data_cluster{end+1,1}=[index_convex1; sub_01]';
           data_cluster{end+1,1}=[index_convex2; sub_02]';
        
        if c_outer == index_convex1;
            sub_outer_01 = k_boundary2(1:i_devide,:);
            sub_outer_02 = k_boundary2(i_devide+1:end,:); 
        elseif c_outer == index_convex2;
            sub_outer_01 = k_boundary2(1:i_devide-1,:);
            sub_outer_02 = k_boundary2(i_devide:end,:); 
        end
        data_outer{end,1}=[index_convex1; sub_outer_01]';%data_outer{end+1,1}=[index_convex1; sub_outer_01]';
        data_outer{end+1,1}=[index_convex2; sub_outer_02]';
            
       end
    end 


    % re-order
    %{
    for j=size(data_cluster,1):-1:2
      aa = ismember(data_cluster{j,1}, data_cluster{j-1,1});
      bb = ismember(data_cluster{j-1,1}, data_cluster{j,1});
      
      if nnz(aa(2:end)) == length(aa(2:end))
          data_cluster(j,:) = [];
      elseif nnz(aa(2:end)) == length(bb(2:end))
          data_cluster(j-1,:) = [];
      end
    end
     %}     
      
end



% Helper function to determine if extreme points are involved
function flag_extreme = getFlagExtreme(index_sub)
    if index_sub(2) == 1 && index_sub(end) ~= 2
        flag_extreme = 1;
    elseif index_sub(2) ~= 1 && index_sub(end) == 2
        flag_extreme = 2;
    elseif index_sub(2) == 1 && index_sub(end) == 2
        flag_extreme = 13;
    else
        flag_extreme = 0;
    end
end        

% Helper function to adjust for extreme points
function b_test2 = adjustExtremePoints(b_test2, flag_extreme, n, size_z)
    if n == 2 && flag_extreme == 1
        b_test2(1, :) = [];
    elseif n == size_z - 1 && flag_extreme == 2
        b_test2(end, :) = [];
    elseif flag_extreme == 13
        b_test2([1, end], :) = [];
    end
end
