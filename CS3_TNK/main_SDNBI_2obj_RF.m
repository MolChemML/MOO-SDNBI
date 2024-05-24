%{
===========================================================================
 Title: Sandwich + mNBI (SDNBI) Algorithm for Multiobjective Optimisation
  main_SDNBI_2obj.m: Main algorithm for bi-objective optimsation applicable
                     for NLP and MINLP problems. The algorithm is mainly
                     based on Lee et al. (2021) that combine sandwich
                     algorithm proposed by Renne et al.(2010) and Normal
                     Boundary Intersection propsed by Das and Dennis (1998)

  Methodology and Code developed by:
   Ye Seol Lauren Lee (ye.lee17@imperial.ac.uk), Imperial College London
   Edward J Graham (edward.graham10@imperial.ac.uk)
               
  Code Log:
    Last update on 14th December 2020
===========================================================================
%}
%% ------------------------------------------------------------------------
%                          0. Initialisation 
%  ------------------------------------------------------------------------
clc; clear;
hold off

% load true Pareto front
% if doens't exist, comment out
load('PPoints_CS3_oNBI.mat'); load('PPoints_CS3_oNBI.mat'); 


tMainStart      = tic;    % Log cpu time (MATLAB simulation)
dim             = 2;      % number of objectives 
n_dim           = 2;      % number of decision variables
n_sobol         = 20;     % number of starting points (generated by sobol' sequnce)
dTol            = 0.005;  % Convergence critera for the algorithm (error between inner and outer)
sTol            = 0.04;   % Convergence criteria for the algorithm (spreadness) 
iter            = 1 ;     % Iteration counter
n_sub           = 1;      % Initialisze number of subproblmes

PPoints         = [];     % store data of Pareto points
facets          = [];     
dontRunIndex    = {};
dontRunPoints   = {};
changeBetaPoints= {};
discard_PPoints = [];

%% -----------------------------------------------------------------------
%             1. Find All Extreme Points (Anchor Points) ZA
%  -----------------------------------------------------------------------
%{
 wMat            = eye(dim); % weights to generate extreme point 
 %wMat(end+1,:)   = ones(1,dim)/norm(ones(1,dim));
 for i = 1:size(wMat,1)
     wMat(i,:) = wMat(i,:)/sum(wMat(i,:));
 end 
   
 mat_subPrb={};

 for i = 1:size(wMat,1)
   
     [objsNew, n_groupsNew, properties, tStore, tElapsed]  = ...
                  f_MtoG_WS_2obj(wMat(i,:),n_sobol);
      

     PPoints(end+1,:) = [ wMat(i,:), objsNew, n_groupsNew, properties, zeros([1,4*dim]), tStore, tElapsed,10^9,1,0];
     mat_subPrb{end+1,1}=[1];
     disp(['Objectives = ' num2str(objsNew)])
   
     filename1 = sprintf('NBI_2obj_cs3_anchor.mat');
     save(filename1);

 end
 tMainElap = toc(tMainStart);
 time_indvidual = cputime;
 filename1 = sprintf('NBI_2obj_cs3_anchor.mat');
 save(filename1);
%}
 load('NBI_2obj_cs3_anchor.mat'); 
 %PPoints(1,end+1:end+3)=[0 0 0]; % to save time I didn't run anchor points

%% -----------------------------------------------------------------------
%        2. Set Parameter and Initialise Variable for SDNBI
%  -----------------------------------------------------------------------
wMat        = PPoints(:,1:dim);
z           = PPoints(:,dim+1:2*dim);
A           = -wMat;
b           = -sum(wMat.*z,2);

dMax        = inf;     % Stopping criteria Initialisation
dHistory    = [];      % To store toleracne history
hvHistory   = [];      % Quality (Hypervolume) history
flag_beta   = 1;
PPointsHistory    = PPoints; % Store anchor points
idx         = 1;       % number of subregions to explore  

% Characterise the point to devide regions in to subregions based on: 
% -3: convex and -9: concave (nonconvex) region
% Initial assumption: Pareto front is defined as convex.
data_subPrbs{1,1} =[-3,1,2]; % Store information of anchor points
data_subOuter{1,1}=[-3,1,2]; % Store information of anchor points
cnt         = 0;  % cluster counter
w_store     = [1 0; 1 0];


%load('CS3_SDNBI_region_210217-findextreme-iter5.mat');

iter=iter+1;
dTol = 0.002;
sTol = 0.035;
%% -----------------------------------------------------------------------
%                  3. Solve Bi-objective Problem by SDNBI
% ------------------------------------------------------------------------
while dMax > dTol
      idx=sum(~cellfun(@isempty,data_subPrbs),2); cnt=1; n1=1;
  for n=n1:idx(end,1)   %length(unique(round(PPoints(:,end)-0.1,0))) % iterations for subproblems (subregions)
      %% ================================================================
      %  3-1) Inner approximation and error measurement
      % =================================================================
      % For the Pareto points, redefine its properties (if it is in
      % convex region: -3, nonconvex region:-9) then construct the convex
      % hull (inner approximation) for subregions 
      % =================================================================
      % Refine the current set based on the covexity assumption imposed on 
      % for convex or concave hull generation
      index_dontadd = 0;
      index_sub = data_subPrbs{iter-1,n}; ifconvex = index_sub(1); flag_beta=1; flag_b=1;%[index_sub, col] = find(PPoints(:,end)==n | round(PPoints(:,end),0)==n | round(PPoints(:,end)-0.1,0)==n);
      index_outer = data_subOuter{iter-1,n};
      [index_outer_corrected] = f_outer_test(PPoints,index_sub,abs(A),abs(b));
    
      z = PPoints(index_sub(2:end), dim+1:2*dim) ; 
      z0 = PPoints(:,dim+1:2*dim); PPoints_temp = PPoints(index_sub(2:end), :);
      dummyMat = generate_dummy(z,dim,ifconvex);
      z_dummy = [z; dummyMat];

      % Generate convex hull (inner approximation)
      [facetIndices, facetNormals] = convexHull_R1(z_dummy,dim); 
      FacetID =  convhulln(z_dummy,{'Qt','n'});
      direction = 1; 
      if ifconvex==-9; direction = -1; end 
      facetNormals = direction * facetNormals;
 
      
      % Plot Parto points with dummy points
      hold off
      k=FacetID;
      xlim([0,1.2]);
      ylim([0,1.2]);
      plot(z_dummy(:,1),z_dummy(:,2),'o','MarkerEdgeColor',rgb('Salmon'))
      hold on
      
      xlim([0,1.2]);
      ylim([0,1.2]);
      plot(PPoints(:,dim+1),PPoints(:,dim+2),'o','MarkerEdgeColor','r')
      plot(PPoints_CS3_oNBI(:,dim+1),PPoints_CS3_oNBI(:,dim+2),'o','MarkerEdgeColor',rgb('LightGray'),'MarkerSize',1.5);
      %plot(PF_CS3_oNBI(:,3),PF_CS3_oNBI(:,4),'o','MarkerEdgeColor',rgb('DodgerBlue'),'MarkerSize',3)
      plot(z_dummy(k,1),z_dummy(k,2),'--','Color','b') %[0, 0.75, 0.75]
 
      % =================================================================
      % Remove facets if: 1) the facet is generated by only a set of dummy
      % points; 2) the normal vector of the facet is outward.
      % =================================================================
      for i = size(facetIndices,1):-1:1 
          facInd = facetIndices{i};
          count = 0;
          for j = 1:length(facInd)
              if facInd(j) <= size(z,1)  
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
      
      % =================================================================
      % Calculate error (btw inner and outer approximations) and crowding 
      % distances (optional) for each facets. The calculation of the error 
      % is based on linear programming (SIMPEX) which identify a vertex that
      % maximise euclidean distance between outer and inner appox.
      % =================================================================
      delta = zeros(size(facetIndices,1),1);
      %z_sort = PPoints(:,3:4); %sortrows(PPoints(:,3:4),1,'ascend'); 
      %space2 = compute_crowding_distances(z_sort');
      space= zeros(size(facetIndices,1),1);
      tInner=0;
      
      for i = 1:size(facetIndices,1) % calculate error for each facet
          flag_dontrun = 0;
          for j = 1:size(dontRunPoints,1)
              PointIndice = facetIndices{i};
              oPtIndice = index_sub(PointIndice+1);
              cell2 = cell2mat(dontRunPoints);
              c_pp=ismember( oPtIndice, dontRunPoints{j} );
            
              if prod(c_pp) == 1;
                 flag_dontrun = 1;
                 break
              end
          end
          for j = 1:size(dontRunIndex,1) % dont run flagged facets
              if facetNormals(i,:) == dontRunIndex{j}
                 flag_dontrun = 1;
                 break
              end
          end
          if flag_dontrun == 1
             delta(i) = 0;
             space(i) = 0;
             continue
          end
          
          iPoints = z_dummy(facetIndices{i},:);
          plot(iPoints(:,1),iPoints(:,2),'o','MarkerEdgeColor','black','LineWidth',2)
        
          %iPoints = z(facetIndices{i},:);
          w = facetNormals(i,:); w_norm=w/sum(w);
          %f = direction * w_norm; 
          f = direction * w;
          options = optimset('linprog');
          options.Display = 'off';
        
          lb = zeros(1,dim);
          lb = min(z);
          ub = max(z);
          index_for_error = index_outer_corrected(1,2:end);
          
          tInner_tic = cputime; % measure the CPU time for error approximation
          zBar = linprog(f,direction*A(index_for_error,:),direction*b(index_for_error,:),[],[],lb,ub,options);
          tInner_mid = (cputime -tInner_tic);
          tInner = tInner + tInner_mid; % measure the CPU time for error approximation;
          
          plot(zBar(1,:),zBar(2,:),'x','MarkerEdgeColor','cyan','LineWidth',2,'MarkerSize',8);
          wT_zBar = w * zBar; 
          ib = w*iPoints(1,:)';
          zUB = max(z0);
          zLB = min(z0);
          epsilon_tolerance = zUB - zLB;

          %delta(i) = abs((ib - wT_zBar)); % error of OPS-IPS
          delta(i) =  abs((ib - wT_zBar)/(w*epsilon_tolerance')); % epsilon-dominance.
          space(i) =  sum(abs(iPoints(1,:)-iPoints(2,:))); %compute_crowding_distances_v2(iPoints');
          
          draw=index_for_error;
          for k=1:length(draw)
              xlim([0, 1.2]); % Fixed x-axis range
              ylim([0, 1.2]); % Fixed y-axis range
              fplot(@(x) -A(draw(k),1)/A(draw(k),2)*x+ (b(draw(k)))/A(draw(k),2),[0,1],'Color','black');
          end
          %fplot(@(x4) -A(4,1)/A(4,2)*x4+ (b(4)/A(4,2)),[0,1],'Color','black');   
      end

      % =================================================================
      % Choose the facet that has maximum error, if the max. error is below
      % the criteria, then space is used.
      % !!! Need to improve this part to have well distributed Pareto front
      % =================================================================    
      if isempty(facetIndices)
         dMax =0;
      else
         [dMax, index] = max(delta);
         sMax  = max(space);
         disp(dMax)
      end
      
      % for the subregion that satisfy the tolerance criteria (inner and
      % outer distance, dMax, and distribution, sMAX), opt out from the
      % consideration
      if (dMax < dTol || sMax < sTol)
         if isempty(facetIndices)
            sMax =0;
         else
            [sMax, index] = max(space); 
             %sMax = 0; % use only dMax as stopping criteria
         end
         if sMax < sTol
            dMax= dTol+1; idx_test=sum(~cellfun(@isempty,data_subPrbs),2); 
            if  size(idx_test,1) >= iter
                continue
            elseif (n == idx_test(iter-1) && size(idx_test,1) < iter)
                    dMax = 0;
                    break;
            else       
                dMax= dTol+0.1;
                continue 
            end
         end
         dMax= dTol+0.1;
      end
      
      % =================================================================
      % C: I think we don't need this part anymore. This was to change
      %    beta manually to get new Pareto points
      % =================================================================
      %{
      facInd = facetIndices{index};
      for j=1:size(changeBetaPoints,1)
          c_beta = ismember(index_sub(facInd+1),changeBetaPoints{j,1});
          if prod(c_beta) == 1
             flag_beta = changeBetaPoints{j,2};
             break
          end
      end
      ====================================================================
      %}
      
      facInd = facetIndices{index};
      params = []; % Plot the facet we are considering at this point.
      tt=[];
      params = PPoints(index_sub(facInd+1),:); % store the objective function that counstruct the facet
      plot( params(:,dim+1), params(:,2*dim),'o','MarkerEdgeColor','black','MarkerFaceColor',rgb('Gold'),'LineWidth',1,'MarkerSize',6);
        
      w = facetNormals(index,:); 
      w_norm =  w./sum(w);
      for k=1:size(w_norm,2)
          if w_norm(1,k) == 0;
             w_norm(1,k) = 10e-5;
          end
      end
      disp(w_norm)
      disp(w)

      %% ================================================================
      %  3-2) RUN mNBI (subproblem of SDNBI) 
      % =================================================================
      % Set parameters of mNBI. n is determined by normalised w, beta is
      % choosen as middle point of the vertex of the hyper plane
      % In 2D beta = [0.5 0.5]
      % =================================================================
      [f_o, i_min] = min(sortrows(params(:,3:4),'ascend'));
      f_k = sortrows(params(:,3:4),'ascend');
    

      if flag_beta == 1;
         %beta = ones(1,size(f_k,1))/sum(ones(1,size(f_k,1))); % middle point
         beta_1= 0.5; 
         beta = [beta_1 1-beta_1]; % end point
      %elseif flag_beta ==3;
      %   beta_1 = 0.1; beta = [beta_1 1-beta_1];
      %elseif flag_beta ==4;
      %   beta_1=0.9; beta = [beta_1 1-beta_1];
      %else
      %   beta_o = rand([1,size(f_k,1)]);
      %   beta = beta_o / sum(beta_o);
      end
      
      Mat_facet = transpose(f_k);
      Mat_ref = repmat(transpose(f_o), 1, dim);
      phi = Mat_facet - Mat_ref;
    
      weight = beta;
      normal = -w_norm; % n in nbi constraints: tn + phi*beta - (f - f_o)
    
      point_on_CHIM = (phi*beta')'+f_o;
      plot(point_on_CHIM(:,1),point_on_CHIM(:,2),'*','MarkerEdgeColor','m');
      fplot(@(x) normal(2)/normal(1)*(x-point_on_CHIM(1))+ point_on_CHIM(2),...
                [point_on_CHIM(1)-0.1,point_on_CHIM(1)+0.1],'Color',rgb('Magenta'));
      
      % Generate integer cut to avoid the generation of identical solutions
      %(This might cuase "covergence fail' if there is no solutions)
      icut = zeros([size(facInd,2) dim]);
      for n_icut=1:size(icut,1)
          icut(n_icut,:)= z(facInd(n_icut),:);
      end
      icut(end,:)=[];

      % =================================================================
      % Solve mNBI in GAMS environment (This part is dependent of the choice
      % of optimisation platform and the problem, which means you need to
      % customisae f_MtoG_NBI_2obj). 
      % =================================================================
      [objsNew, n_groupsNew,properties,tStore, tElapsed] = ...
        f_MtoG_NBI_2obj(weight, beta, normal, phi, f_o, f_k, icut, dim, n_dim, n_sobol);
    
      PPointsHistory(end+1,:) = [ w, objsNew, n_groupsNew, properties, beta, point_on_CHIM, tStore, tElapsed, dMax, sMax, tInner];
      
 
      % =================================================================
      % Investigate the constraints: phi*beta + tn >= f - f_o
      % if any of these constratints are inactive, backtrack the
      % coressponding beta and screen out the facet, as there will be no
      % new Pareto points. The next normal vector is also achieved by
      % a set of Lagrange multiplier for the constraints 
      % isasctive =1 if all of those constraints are active; otherwise 0
      % =================================================================
      nbi_const = zeros([1 dim]);
      nbi_const = round(properties(1,2:dim+1),5);
        
      isactive = ~any(nbi_const);
      isidentical = 0;
      diffMat = bsxfun(@minus,z,objsNew);
      diff = (sum(diffMat.^2,2)).^0.5;
      if any(diff < 1e-5)
         % Flag not to be chosen anymore
         disp('Pareto Points already found: Go to mNBIn(a)/(b)')
         isidentical = 1;
      end
         
      if (isactive == 0 && isidentical == 0)  % Backcalcuation of the beta that preserve the nbi equation equality for biobjective prbs
         syms beta_new t;
         [sol_beta,sol_t]= solve( phi(1,1)*beta_new + phi(1,2)*(1-beta_new)+t*normal(1)-(objsNew(1)-f_o(1))==0, ...
                 phi(2,1)*beta_new + phi(2,2)*(1-beta_new)+t*normal(2)-(objsNew(2)-f_o(2))==0 );           
  
          z_nopt =double(phi*[sol_beta;1-sol_beta])'+f_o;
          plot(z_nopt(1),z_nopt(2),'x','MarkerEdgeColor','g');
         
          beta_new = double([sol_beta, 1-sol_beta]); % this is to get lagrange multiplier (make the nbi equality constraint active)
          [objsNew_e, n_groupsNew_e,properties_e,tStore_e, tElapsed_e] = ...
               f_MtoG_oNBI_2obj(weight, beta_new, normal, phi, f_o, f_k, icut, dim, n_dim, n_groupsNew);
%           properties(1) = properties_e(1);
%           properties(1,4) = properties_e(1,4);
%           properties(1,5) = properties_e(1,5);
%           properties(1,6) = properties_e(1,6);
%           properties(1,7) = properties_e(1,7);
      end
       
      % store optimal parameters    
      z_star              = objsNew; 
      n_groups_star       = n_groupsNew;
      properties_star     = properties;

        
      %% ================================================================
      %  3-3) Check if the new Pareto Point is identical to the others
      % =================================================================
      % Check if the identified pareto (z_star) point already chosen
      % If the z_start is found to be identical to the exsiting Pareto
      % Points, solve bilevel programming to search new normal vector that
      % results in new Pareto Points. 
      % =================================================================     
      diffMat = bsxfun(@minus,z0,z_star);
      diff = (sum(diffMat.^2,2)).^0.5;
      if any(diff < 1e-5)
         % Flag not to be chosen anymore
         disp('Pareto Points already found: Solve bilevel')
        %  [objsNew, n_groupsNew, properties, normalNew, tWall, tCPU] ...
        %= f_MtoG_NBI_2obj_normal(weight, beta, normal, phi, f_o, icut,dim,n_dim, n_sobol,iPoints);
        % change beta to find the new point for same facet
        if beta(1) == 0.9 || beta(1) == 0.1
           dontRunIndex{end+1,1} = w; %need to change this beta !!!!!
           dontRunPoints{end+1,1}= index_sub(facInd+1);
           discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, i_sub, tInner];
  
           continue
        end
         
        % Determine how we vary nomral vectore (clockwise or anti-clockwise for biobjective)      
         [f1_min,id3]=min(params(:,3));
         [f2_min,id4]=max(params(:,3));
        if round(objsNew(1),5)== round(f1_min(1),5)
            flag_b = 3;  % beta = 0.5~1 has no solution : need to rotate normal vector anti-clockwise
            log_objs(1,:) = objsNew; what_to_exp=params(id4,3:4);
        else
            flag_b = 4;  % beta = 0~0.5 has no solution : need to rotate normal vector clockwise
            log_objs(1,:) = objsNew; what_to_exp=params(id3,3:4);
        end      

        
        %-----------------------------------
        %  Solve NBIn: simple version
        %-----------------------------------
        
        %[objsNew, n_groupsNew, properties, normalNew, tWall, tCPU] ...
        %    = f_MtoG_NBI_2obj_normal(weight, beta, normal, phi, f_o, icut,dim,n_dim, n_sobol,iPoints,flag_b);

        [objsNew, n_groupsNew, properties_e, tStore_e, tElapsed_e] ...
            = f_MtoG_mNBIn_2obj(weight, beta, normal, phi, f_o, f_k, icut,dim,n_dim, n_sobol, flag_b);
        

      
       %fplot(@(x) normalNew(2)/normalNew(1)*(x-point_on_CHIM(1))+ point_on_CHIM(2),...
       %         [point_on_CHIM(1)-0.1,point_on_CHIM(1)+0.1],'Color',rgb('Salmon'));
            
        syms beta_new t;
       [sol_beta,sol_t]= solve( phi(1,1)*beta_new + phi(1,2)*(1-beta_new)+t*normal(1)-(objsNew(1)-f_o(1))==0, ...
                    phi(2,1)*beta_new + phi(2,2)*(1-beta_new)+t*normal(2)-(objsNew(2)-f_o(2))==0 );           

        z_nopt =double(phi*[sol_beta;1-sol_beta])'+f_o;
        plot(z_nopt(1),z_nopt(2),'x','MarkerEdgeColor','g');
        fplot(@(x) normal(2)/normal(1)*(x-z_nopt(1))+ z_nopt(2),...
                [z_nopt(1)-0.1,z_nopt(1)+0.1],'Color','g'); 
        plot(objsNew(:,1),objsNew(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r')
        
        beta_new = double([sol_beta, 1-sol_beta]); % this is to get lagrange multiplier (make the nbi equality constraint active)
   

        PPointsHistory(end+1,:) = [ 99,99 , objsNew, n_groupsNew, properties_e, beta_new, point_on_CHIM, tStore_e, tElapsed_e, dMax, sMax, 0]; 
        
                
        %-----------------------------------
        %  Solve NBIn: Bilevel formulation
        %-----------------------------------
        %{
       syms n_new t;
       [sol_n,sol_t]= solve( phi(1,1)*beta(1) + phi(1,2)*beta(2)+t*n_new-(what_to_exp(1)-f_o(1))==0, ...
                    phi(2,1)*beta(1) + phi(2,2)*(beta(2))+t*(-1-n_new)-(what_to_exp(2)-f_o(2))==0 );           
        n_new = double([sol_n, -1-sol_n]);
        fplot(@(x) n_new(2)/n_new(1)*(x-point_on_CHIM(1))+ point_on_CHIM(2),...
                [point_on_CHIM(1)-0.1,point_on_CHIM(1)+0.1],'Color','g'); 

            
        n_groups_ini = params(:,2*dim+1:2*dim+n_dim);
        [objsNew, n_groupsNew, properties, normalNew, tWall, tCPU] ...
            = f_MtoG_NBI_2obj_normal(weight, beta, normal, phi, f_o, icut,dim,n_dim, n_sobol,params(:,3:4),n_groups_ini,flag_b);
        
      
       fplot(@(x) normalNew(2)/normalNew(1)*(x-point_on_CHIM(1))+ point_on_CHIM(2),...
                [point_on_CHIM(1)-0.1,point_on_CHIM(1)+0.1],'Color',rgb('Salmon'));
            
        syms beta_new t;
       [sol_beta,sol_t]= solve( phi(1,1)*beta_new + phi(1,2)*(1-beta_new)+t*normal(1)-(objsNew(1)-f_o(1))==0, ...
                    phi(2,1)*beta_new + phi(2,2)*(1-beta_new)+t*normal(2)-(objsNew(2)-f_o(2))==0 );           

        z_nopt =double(phi*[sol_beta;1-sol_beta])'+f_o;
        plot(z_nopt(1),z_nopt(2),'x','MarkerEdgeColor','g');
        fplot(@(x) normal(2)/normal(1)*(x-z_nopt(1))+ z_nopt(2),...
                [z_nopt(1)-0.1,z_nopt(1)+0.1],'Color','g'); 
        plot(objsNew(:,1),objsNew(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r')
        
        beta_new = double([sol_beta, 1-sol_beta]); % this is to get lagrange multiplier (make the nbi equality constraint active)
   
        [objsNew_e, n_groupsNew_e,properties_e,tStore_e, tElapsed_e] = ...
            f_MtoG_oNBI_2obj(weight, beta_new, normal, phi, f_o, f_k, icut, dim, n_dim, round(n_groupsNew,7));

        PPointsHistory(end+1,:) = [ 99,99 , objsNew, n_groupsNew, properties, beta_new, point_on_CHIM, tStore, tElapsed, dMax, sMax, 0]; 
        %}
            
        properties(1) = properties_e(1);
        properties(1,4) = properties_e(1,4);
        properties(1,5) = properties_e(1,5);
        properties(1,6) = properties_e(1,6);
        properties(1,7) = properties_e(1,7);
       % END of NBIn
       
        z_star              = objsNew; % store optimal parameters
        n_groups_star       = n_groupsNew;
        properties_star     = properties;
        flag_b=999; % 999: log dontRunPoints  
        log_objs(2,:) = z_star;
        
        if flag_b ==3 || flag_b==4;
        changeBetaPoints{end+1,1}=index_sub(facInd+1);
        changeBetaPoints{end,2}= flag_b;
        
        discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, i_sub, tInner];
        data_subPrbs{iter,cnt}=index_sub; cnt=cnt+1;
        continue
        end
      end % END of BILEVEL
      flag_beta = 1;
  
      %% ================================================================
      %  3-5) Check the quality set of Praeto Points and Updat IPS
      % =================================================================
      % Check if the ident
      % =================================================================       
      iPoints = z_dummy(facetIndices{index},:);
      ib = w * iPoints(1,:)';
      [pFront,idxs] = paretoFront([z;z_star]*-1);  
      pFront = -pFront;
      
      diffMat = bsxfun(@minus,z0,z_star);  %diffMat = bsxfun(@minus,z,z_star);
      diff = (sum(diffMat.^2,2)).^0.5;
      if any(diff < 1e-5)
         disp('Pareto Points already found: Don Run This Facet')
         % Flag not to be chosen anymore
         dontRunIndex{end+1,1} = w; %need to change this beta !!!!!
         index_dontadd = 1;
         dontRunPoints{end+1,1}= index_sub(facInd+1);
         discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, i_sub, tInner];
      end
        
      % Improve this part. (not in use untill 100 Pareto points are
      % generated)
      if size(PPoints,1) >= 100
         if idxs(end) ~= size(z,1) + 1
            disp('new PPoint not Pareto optimal, flag this facet not to run again')
            dontRunIndex{end+1,1} = w;
            discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, i_sub, tInner];
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
            discard_PPoints(end+1,:)=[iter, w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, i_sub, tInner];
            continue
         end
      end

      % =================================================================
      % Update inner approximation (IPS) by adding new point to the set of 
      % Pareto points if new points is not dominated by others. 
      % =================================================================  
      [M,i_extreme]= min(z(:,1));
      [M,I]= min(z(:,2));
      i_extreme(end+1)=I;
      z = [z;z_star];
      plot(z_star(:,1),z_star(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r')
      i_sub = round(PPoints_temp(facInd(1),end)-0.1,0);
   
      [ndf, df]=non_dominated_front([PPoints(:,3:4);z_star]');
      if ~isempty(df)
         disp('new point found to be dominated')
      end
      
      if index_dontadd == 0;
        PPoints(end+1,:) = [ w, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed, dMax, sMax, tInner];
      

      %% ================================================================
      %  3-6) Update outer approx (OPS) by adding inequality wTz <= b
      % =================================================================
      % check if new point is unsupproted (nonconvex), if then, wTz>=b &
      % Update the index for subregions and cluster or devide the regions 
      % to subregions based on the its properties (convex vs. nonconvex)
      % =================================================================      
      multiplier= -1;  
    %z2(end+1,:)=[z2(1,1),z2(1,2)+3]; z2(end+1,:)=[z2(end-1,1)+3,z2(end-1,2)]; z2(end+1,:)=[z2(end,1),z2(end-1,2)];
    %if properties(1) <0 
    %    multiplier =1;
    %end
  
      w_correct =  PPoints(end,12:13)/sum(PPoints(end,12:13)); % from Lagrangian multiplier
      w_store(end+1,:)=w_correct;
      A(end+1,:) = multiplier*w_correct; %multiplier*PPoints(end,1:dim);
      b(end+1,:) = multiplier*sum(w_correct.*PPoints(end,dim+1:2*dim),2);
    
      fplot(@(x) -A(end,1)/A(end,2)*x+ (b(end))/A(end,2),[0,1],'Color',rgb('OrangeRed'));
      %normal vector
      fplot(@(x) w_correct(2)/w_correct(1)*(x-z_star(1))+ z_star(2),...
                [z_star(1)-0.1,z_star(1)+0.1],'Color',rgb('OrangeRed'));
            
      c_outer = -3;
      if properties(1) < 0
         c_outer = -9;
      end
     [data_cluster, b_test2, data_outer] = f_Cluster(index_sub,PPoints,A,b,z_star,dim, c_outer);
     [data_cluster9, b_test9, data_outer9] = clusterParetoPoints(index_sub,PPoints,A,b,z_star,dim, c_outer);
     %[data_cluster5, b_test5, data_outer5] = f_cluster_Pareto_v5(index_sub,PPoints,A,b,z_star,dim, c_outer);
     [data_cluster1, b_test, data_outer1] = f_cluster_Pareto_v3(index_sub,PPoints,A,b,z_star,dim, c_outer);
     [data_cluster2, b_test, data_outer2] = f_cluster_Pareto_v1(iPoints,z,PPoints,b(end,1),multiplier,w_correct,dim,c_outer);
     disp(b_test'); 
     %A(end,:)= multiplier*A(end,:); b(end,1) = multiplier*b(end,:);
    
     % Store data
     for j=1:size(data_cluster,1)
         data_subPrbs{iter,cnt} = data_cluster{j,1};
         data_subOuter{iter,cnt} = data_outer{j,1};
         cnt = cnt + 1;
     end
     end % END of index_dontadd
     % Opt out facet of which no further Pareto points exist.
     % Need to improve for more general case
     if (flag_b==999)
         i_match = [0, 0];
         i_match1 = find( round(log_objs(1,1),5) == round(PPoints(:,dim+1),5) );
         if length(i_match1)>1
            i_match1 = i_match1(find( round(log_objs(1,2),5) == round(PPoints(i_match1,dim+2),5) ));
         end
         i_match2 = find( round(log_objs(2,1),5) == round(PPoints(:,dim+1),5) );
          if length(i_match2)>1
            i_match2 = i_match2(find( round(log_objs(2,2),5) == round(PPoints(i_match2,dim+2),5) ));
          end
         i_match = [i_match1, i_match2];
         dontRunPoints{end+1,1}= i_match; %index_sub(facInd+1);  
         flag_b = 1;
     end
  
     dHistory(end+1) = dMax;
     %PPointsHistory(end+1,:) = PPoints(end,:);
     
     
     filename = sprintf('CS3_SDNBI_region_UROP-Tol=0.05_iter.mat',iter);
     save(filename);
     tMainElapsed=toc(tMainStart);
     
 end
  save('CS3_SDNBI_region_UROP.mat');
  iter=iter+1; 
end

%% -----------------------------------------------------------------------
%                 4. Evaluate final results and Save it
% ------------------------------------------------------------------------
PPoints_CS3_SDNBI = PPoints;
z_final = PPoints(:,3:4);
f_u = [0.0416641269, 0.0416641269]; f_n = [1.0384498374, 1.0384498374];
z_final(:,1)=(f_n(1)-f_u(1))*z_final(:,1)+f_u(1);
z_final(:,2)=(f_n(2)-f_u(2))*z_final(:,2)+f_u(2);
hv_CS3_SDNBI=approximate_hypervolume_ms(PPoints(:,3:4)',[1,1]');

[ndf_index, df_index] = non_dominated_front(CS3_Analysis_mNBI(1:58,2:3)');
PF_CS3_mNBI = CS3_Analysis_mNBI(ndf_index,:);

[ndf_index, df_index] = non_dominated_front(PPoints_CS3_SDNBI(:,3:4)');
PF_CS3_SDNBI = PPoints_CS3_SDNBI(ndf_index,:);


[ndf_index, df_index] = non_dominated_front(PPoints(:,3:4)');
PF_CS3_SDNBI = PPoints(ndf_index,:);

crowding_dist=compute_crowding_distances(PPoints(:,3:4)');
CD_CS3_SDNBI=crowding_dist;
CD_stdv_CS3_SDNBI=std(CD_CS3_SDNBI(:,dim+1:end)');
CD_avg_CS3_SDNBI =mean(CD_CS3_SDNBI(:,dim+1:end)');
DM_CS3_SDNBI = compute_quality_pareto(PPoints(:,3:4), [1,1], [0,0]);
%save('CS3_SDNBI_region_210217-R5-crowding_distance.mat');

% For plotting ------
hold off;
plot(PPoints_CS3_oNBI(:,dim+1),PPoints_CS3_oNBI(:,dim+2),'o','MarkerEdgeColor',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',6);
hold on;
plot(PF_CS3_oNBI(:,dim+1),PF_CS3_oNBI(:,dim+2),'o','MarkerEdgeColor',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',6);
plot( PF_CS3_SDNBI(:,3), PF_CS3_SDNBI(:,4),'o','MarkerEdgeColor','red','MarkerFaceColor',rgb('Gold'),'LineWidth',1,'MarkerSize',6);
plot( PF_CS3_SD(:,3), PF_CS3_SD(:,4),'x','MarkerEdgeColor','blue','MarkerFaceColor',rgb('blue'),'LineWidth',1,'MarkerSize',6);
plot( PF_CS3_mNBI(:,2), PF_CS3_mNBI(:,3),'^','MarkerEdgeColor','green','MarkerFaceColor',rgb('green'),'LineWidth',1,'MarkerSize',6);


DM_CS3_SDNBI = compute_quality_pareto(PF_CS3_SDNBI(:,3:4), [1,1], [0,0]);
DM_CS3_SD = compute_quality_pareto(PF_CS3_SD(:,3:4), [1,1], [0,0]);
DM_CS3_mNBI = compute_quality_pareto(PF_CS3_mNBI(:,3:4), [1,1], [0,0]);
DM_CS3_oNBI = compute_quality_pareto(PF_CS3_oNBI(:,3:4), [1,1], [0,0]);

hv_CS3_SDNBI=approximate_hypervolume_ms(PF_CS3_SDNBI(:,3:4)',[1.1,1.1]');
hv_CS3_SD=approximate_hypervolume_ms(PF_CS3_SD(:,3:4)',[1.1,1.1]');
hv_CS3_mNBI=approximate_hypervolume_ms(PF_CS3_mNBI(:,3:4)',[1.1,1.1]');
hv_CS3_oNBI=approximate_hypervolume_ms(PF_CS3_oNBI(:,3:4)',[1.1,1.1]');


