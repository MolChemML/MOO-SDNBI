%%---------------- HIGH DEMENTIONAL SANDWICH ALGORITHM -----------------%%
%{
THIS CODE IS DEVELOPED BASED ON EDWARD GRAHAM MALAB CODE(27-MAR-2018)
The algorithm is taken from Rennen(2009) - 

%}
clc; clear;
%%-------------- GPROMS SETUP FOR MINLP OPTIMISATION --------------------%%
tMainStart      = tic;
dim             = 2;      % number of objectives 
n_sobol         = 20;
dTol            = 0.00001;  % Convergence critera for the algorithm
iter            = 1 ;
%% Initialise Array

PPoints         = []; % [weights(:,1:dim), objs(:,dim+1:2*dim), AADs(:,3*dim + 1:4*dim), parameters (:,4*dim + 1:end)]
PPoints_History = [];
facets          = []; % [indices (1:dim), normals(dim+1:2*dim), boundaryIndex(end-1), d(end)]
dontRunIndex    = {};
discard_PPoints = [];
%% 1. Find All Extreme Points (Anchor Points) ZA
wMat            = eye(dim); % weights to generate extreme point 
wMat(end+1,:)   = ones(1,dim)/norm(ones(1,dim));

for i = 1:size(wMat,1)
    %wMat(i,:) = wMat(i,:)/norm(wMat(i,:));
    wMat(i,:) = wMat(i,:)/sum(wMat(i,:));
end 

SuccessCriteria = string(['Optimal Solution Found']);

% To reduce calculation inactivate this part(find anchor point) and
%dirrectly loaded from previously simulated file.

%{
for i = 1:size(wMat,1)
   
    [objsNew, n_groupsNew, properties, tStore, tElapsed]  = ...
                 f_MtoG_WS_2obj(wMat(i,:),n_sobol);
      
    
    PPoints(end+1,:) = [ wMat(i,:), objsNew, n_groupsNew, properties, tStore, tElapsed, 999,0];
    PPoints_History(end+1,:) = [ wMat(i,:), objsNew, n_groupsNew, properties, tStore, tElapsed, 999,0];
    disp(['Objectives = ' num2str(objsNew)])
   
    filename1 = sprintf('WS_2obj_anchor.mat');
    save(filename1);

end
tMainElap = toc(tMainStart);
time_indvidual = cputime;
filename1 = sprintf('WS_2obj_anchor.mat');
save(filename1);
%}
load('NBI_2obj_cs2_anchor.mat'); PPoints(:,11:19)=[]; dTol=1e-5;


%% 2. Set OPS (OuterApproximated Points)
wMat    = PPoints(:,1:dim);
z       = PPoints(:,dim+1:2*dim);
A       = -wMat;
b       = -sum(wMat.*z,2);

%% 3. Calculate of Each Facets
dMax = inf;  % Criteria Initialisation
dHistory = [];
PPointsHistory = PPoints;

flag_beta =1;

iter=1;
while dMax > dTol
    z = PPoints(:,dim+1:2*dim); ifconvex=-3;
    dummyMat = generate_dummy(z,dim,ifconvex); 
    z_dummy = [z;dummyMat];
    %z_test = z(1:2,1:2);
    %z_test(end+1,:)=[3.34E-01	4.08E-01];
    %dim_test = 2;
    [facetIndices, facetNormals] = convexHull_R1(z_dummy,dim); % IPS defined by indices and normals
    FacetID =  convhulln(z_dummy,{'Qt','n'});
    
      hold off
      plot(PPoints(:,dim+1),PPoints(:,dim+2),'o','MarkerEdgeColor','r')
      hold on
    
    for i = size(facetIndices,1):-1:1 % LL:Determine whether the set is consisted with all dummypoint (>z_ub)
        facInd = facetIndices{i};
        count = 0;
        for j = 1:length(facInd)
            if facInd(j) <= size(z,1)  % 
                count = count + 1;
            end
        end
        if count < 2 % if only dummy points
            facetIndices(i,:) = [];
            facetNormals(i,:) = [];
        end
    end                     
    for i = size(facetNormals,1):-1:1 % remove outward facing facets
        if any(facetNormals(i,:) < 0)
            facetIndices(i,:) = [];
            facetNormals(i,:) = [];
        end
    end
    delta = zeros(size(facetIndices,1),1); tInner=0;
    for i = 1:size(facetIndices,1) % calculate error for each facet
        flag = 0;
        for j = 1:size(dontRunIndex,1) % dont run flagged facets
            if facetNormals(i,:) == dontRunIndex{j}
                flag = 1;
                break
            end
        end
        if flag == 1
            delta(i) = 0;
            continue
        end
        iPoints = z_dummy(facetIndices{i},:);
        w = facetNormals(i,:);
        f = w;
        options = optimset('linprog');
        options.Display = 'off';
        lb = zeros(1,dim);
        tInner_tic=cputime;
        zBar = linprog(f,A,b,[],[],lb,[],options);
        tInner_mid=cputime - tInner_tic;
        tInner=tInner+tInner_mid;
        wT_zBar = w * zBar;
        ib = w*iPoints(1,:)';
        zUB = max(z);
        zLB = min(z);
        epsilon_tolerance = zUB - zLB;
   
        delta(i) = abs((ib - wT_zBar)/(w*epsilon_tolerance'));
      
             
    end

%% 4. Select A Facet with the Largest Error
    [dMax, index] = max(delta);
    disp(dMax)
    if dMax < dTol
        break
    end
 %% 5. Determine New Point Z* By Solving z* argmin{wTz]
    facInd = facetIndices{index};
    params = []; %Nunmber of groups
    
    
    w = facetNormals(index,:);  %% Do I Need to change Facet Normal???
    w_norm =  w./sum(w);
    for k=1:size(w_norm,2)
       if w_norm(1,k) == 0;
           w_norm(1,k) = 10e-5;
       end
    end
    disp(w_norm)
    disp(w)
    
    [objsNew, n_groupsNew,properties,tStore, tElapsed] = ...
        f_MtoG_WS_2obj(w_norm,n_sobol);
    
    z_star              = objsNew; % store optimal parameters
    n_groups_star       = n_groupsNew;
    properties_star     = properties;
    %modelResults_star   = modelResultsBigCell{minIndex};
        
   PPoints_History(end+1,:)= [w, z_star, n_groups_star, properties_star, tStore, tElapsed, dMax, tInner];
    %% RUN CHECKS ON NEW PARETO POINT BEFORE ADDING
    % Check if pareto point already chosen
    
    diffMat = bsxfun(@minus,z,z_star);
    diff = (sum(diffMat.^2,2)).^0.5;
    if any(diff < 1e-7)
         % Flag not to be chosen anymore
        disp('Pareto Points already found')
        %dontRunIndex{end+1,1} = w; %need to change this beta !!!!!
       % if flag_beta < size(beta_store,1)
       %   flag_beta = flag_beta + 1;
       % else
            dontRunIndex{end+1,1} = w;
         %   flag_beta =1;
        %end
        
        discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, tStore, tElapsed,dMax, tInner];
        continue
    end 
    flag_beta = 1;
%% 6. If wTz* = b, Set Error of This Facet to Zero (dont run again)
    iPoints = z_dummy(facetIndices{index},:);
    ib = w * iPoints(1,:)';
    [pFront,idxs] = paretoFront([z;z_star]*-1);  
     pFront = -pFront;
     
     %PPointTemp = PPoints(index,:);
     %PPoints=[];
     %PPoints = PPointTemp;
   
     

    if size(PPoints,1) >= 60
    if idxs(end) ~= size(z,1) + 1
        disp('new PPoint not Pareto optimal, flag this facet not to run again')
        dontRunIndex{end+1,1} = w;
        discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, tStore, tElapsed, dMax, tInner];
        continue
    end
    if length(idxs) ~= size(z,1) + 1
        disp('previous pareto point(s) found not to be optimal, removing these points and constraints')    
        idxs = idxs(1:end-1);
        for k=1:size(PPoints,1)
            if k ~= idxs
                discard_PPoints(end+1,:) = [iter, PPoints(k,:)];
            end
        end
        PPoints = PPoints(idxs,:);
        
        
    end
    if abs(z_star * w' - ib) < 1e-20
        disp('wT*zbar = b, new point found to lie on current facet, flag this facet not to run again')
        dontRunIndex{end+1,1} = w;
        discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, tStore, tElapsed, dMax, tInner];
        continue
    end
    end
    % Add new point to Pareto Front      
                                   
    z = [z;z_star];
    PPoints(end+1,:) = [w, z_star, n_groups_star, properties_star, tStore, tElapsed, dMax, tInner];
    
    % Check if all points are optimal by checking if the weighted sum
    % objective for each point is minimum when giving the same weighting to
    % all other pareto points
    weights = PPoints(:,1:dim);
    points = PPoints(:,dim+1:2*dim);
    
    if size(points,1) > 70
       i = size(weights,1);
           logic1 = [weights(i,:) == weights];
           logic2 = [sum(logic1,2) == dim];
           logic2(i,:)= false;
           if sum(logic2) > 0;
              [maxVal,index] = max(logic2);
              objValue = sum((points.*weights),2);
              if objValue(i,:) > objValue(index,:);
                  discard_PPoints(end+1,:)=[iter,PPoints(i,:)];
                  PPoints(i,:) = [];
                 
              elseif objValue(i,:) < objValue(index,:);
                  discard_PPoints(end+1,:)=[iter,PPoints(index,:)];
                  PPoints(index,:)=[];
               end
           end
      end
           
    %{          
    if size(PPoints,1) > 70
    for i = size(weights,1):-1:1
        objFun = sum(bsxfun(@times,weights(i,:),points),2);
        [v,index] = min(objFun);
       if i ~= index
            disp(fprintf('Pareto point no. %d removed because the same weight gives a better fObj for PPoint %d',i,index))
            PPoints(i,:) = [];
                       dontRunIndex{end+1,1} = w;
         
        end
    end
    end
      %}
%  ASK TO ED
    %% 7. Update IPS by replacing it with conv{z*, IPS} (do this above)
    
    %% 8. Update OPS by adding inequality wTz <= b
    A = -PPoints(:,1:dim);
    b = -sum(PPoints(:,1:dim).*PPoints(:,dim+1:2*dim),2);
    dHistory(end+1) = dMax;
    PPointsHistory(end+1,:) = PPoints(end,:);
   
     

     iter=iter+1;
      tMainElapsed=toc(tMainStart);

end

save('CS2_SD_final_run.mat');

%% -----------------------------------------------------------------------
%                 4. Evaluate final results and Save it
% ------------------------------------------------------------------------
% Run main_Analysis_xx.m