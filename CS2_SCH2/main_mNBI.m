% MultiObjective Optimisation - mNBI method
% Reference points BetaPhi are predefined.
% Updated: 21-AUG-2023 by Lauren LEE
hold off
clc;clear;
%% Optimisation Inupts
PPoints         = []; % Store solutions
discard_PPoints = [];
dim = 2;              % number of objectives
n_dim = 2;            % number of integer control variables
n_sobol =20;          % multi-start
weight = zeros(1,dim);
tCPUinitial = cputime;
tWallinitial = tic;


n_point=10; % initial reference points to start with
n_iter=5;  
%% Main: mNBI method

% -------------------------------------------------------------------------
%  Step 1: Find extreme points
% -------------------------------------------------------------------------
%{
wMat            = eye(dim); % weights to generate extreme point 

for i=1:size(wMat,1)
     [objsNew, n_groupsNew, properties, tStore, tElapsed]  = ...
                 f_MtoG_WS_2obj(wMat(i,:),n_sobol);
     PPoints(end+1,:) = [ wMat(i,:), objsNew, n_groupsNew, properties, zeros([1,4*dim]), tStore, tElapsed,];
end    
%}
 %Instead of running above to get the achor point, here we are loading the
 %saved simulation.
 load('NBI_2obj_cs2_anchor.mat'); 
 PPoints(:,end-2:end) = [];
 
    plot(PPoints(1:dim,1),PPoints(1:dim,2),'o','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',7)
    hold on
% -------------------------------------------------------------------------
%  Step 2: Assign the paramaters
% -------------------------------------------------------------------------

  params=PPoints(1:2,dim+1:2*dim);
  [f_o, i_min] = min(sortrows(params(:,1:2),'ascend'));
    %f_o = zeros([1, dim]);
     
  f_k = sortrows(params(:,1:2),'ascend');
  phi = [f_k(1,1)-f_o(1) f_k(2,1)-f_o(1);...
           f_k(1,2)-f_o(2) f_k(2,2)-f_o(2)];
    
  delta_beta = 1/(n_point-1);
  data_beta=[delta_beta 1-delta_beta];
  for i=2:(n_point-1)
      beta_1= i*delta_beta;
      if beta_1 >= 1
          break;
      end
      data_beta(end+1,:)= [beta_1 1-beta_1];     
  
  end
   
  data_beta = sortrows(data_beta,1,'ascend');  
  flag_ascend = 0;
  normal =-abs((PPoints(2,dim+1:2*dim)-PPoints(1,dim+1:2*dim)))/sum(abs(PPoints(2,dim+1:2*dim)-PPoints(1,dim+1:2*dim)));
  
PPoints_History = PPoints;

%load('CS5_mNBI_210217-3P37.mat');

hold off;

%hold off
%load('CS4_mNBI_210217-1P10-5PXX-nonsobol.mat');
%plot(PPoints(:,3),PPoints(:,4),'o','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',7)
hold on

% -------------------------------------------------------------------------
%  Step 3: Generate Pareto points
% ------------------------------------------------------------------------ 
store_data=[1 0 1; 0 1 1];
i_index = 1;
%load('CS5_mNBI_210217-1P10-3PXX.mat');

for j=1:5
    
    for i=i_index:size(data_beta,1)
        i=i_index; 
        beta = data_beta(i,:); store_data(end+1,:)=[beta,1]; flag_ineq=1;
        disp(beta)
        %weight = beta;  

        %if i==15 || i==16 
        %    i_index= i_index+1;
        %  continue
        %end

        point_on_CHIM = (phi*beta')'+f_o;
        plot(point_on_CHIM(:,1),point_on_CHIM(:,2),'x','MarkerEdgeColor',rgb('DarkBlue'));
        hold on

        icut=[1 1];

        [objsNew, n_groupsNew,properties,tStore, tElapsed] = ...
        f_MtoG_NBI_2obj(weight, beta, normal, phi, f_o, f_k, icut, dim, n_dim, n_sobol);
    
    
        if (round(properties(2),5)~=0) || (round(properties(3),5)~=0)  % Backcalcuation of the beta that preserve the nbi equality
            %t_value = [properties(1,1), properties(1,1)];
            syms beta_new t;
            [sol_beta,sol_t]= solve( phi(1,1)*beta_new + phi(1,2)*(1-beta_new)+t*normal(1)-(objsNew(1)-f_o(1))==0, ...
                                  phi(2,1)*beta_new + phi(2,2)*(1-beta_new)+t*normal(2)-(objsNew(2)-f_o(2))==0 );           
            z_nopt =double(phi*[sol_beta;1-sol_beta])'+f_o;
             plot(z_nopt(1),z_nopt(2),'x','MarkerEdgeColor','g','MarkerSize',7);

             criteria = double(sol_beta);
             if flag_ascend == 0;
                 whereToStart = find(data_beta(:,1)<=criteria) ;
                 i_index = whereToStart(1,1)-1;
                 
             else                
                 whereToStart = find(data_beta(:,1)>=criteria) ;
                 i_index = whereToStart(1,1)-1;
             end
                
                if criteria > beta(1); 
                    store_data(end,:)=[criteria,1-criteria,1]; 
                    store_data(end+1,:)=[beta(1),1-beta(1),0];
                    
                    row = find(round(store_data(1:end-2,1),5)==round(beta(1),5));
                    if ~isempty(row) && round(store_data(row+1,1),5) <= round(criteria,5) && store_data(row+1,3)==0
                        store_data(row:row+1,:)=[];
                    end
                else
                    store_data(end,3)=1; 
                    store_data(end+1,:)=[criteria,1-criteria,0];
                    row = find(round(store_data(1:end-2,1),5)==round(criteria,5));
                    if ~isempty(row) && round(store_data(row-1,1),5) <= round(beta(1),5) && store_data(row,3)==0
                        store_data(row-1:row,:)=[];
                    end
                end
            
             if (i_index <= i && flag_ascend==0) || (i_index >= i && flag_ascend==1)    % for ascend
                 i_index=i;
             end

        end

        z                   = PPoints(:,dim+1:2*dim);
        z_star              = objsNew; 
        n_groups_star       = n_groupsNew;
        properties_star     = properties;  


        diffMat = bsxfun(@minus,z,z_star);
        diff = (sum(diffMat.^2,2)).^0.5;
        if any(diff < 1e-5)
             % Flag not to be chosen anymore
            disp('Pareto Points already found')
            PPoints_History(end+1,:) = [ weight, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed];

            discard_PPoints(end+1,:)=[weight, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed];
            i_index = i_index + 1;
            continue
        end 

        PPoints(end+1,:) = [ weight, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed];
        PPoints_History(end+1,:) = [ weight, z_star, n_groups_star, properties_star, beta, point_on_CHIM, tStore, tElapsed];

        plot(z_star(:,1),z_star(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',7)


        i_index = i_index + 1;
        if i_index > size(data_beta,1)
            break
        end
    end
    
    store_data = sortrows(store_data,1,'ascend');    
    data_beta1=[];
    for i=1:size( store_data,1)-1 %%i=1:size(PPoints,1)-1 
     if store_data(i+1,3) ==1
        bb = 0.5*(store_data(i,1:2) + store_data(i+1,1:2));
        data_beta1(end+1,:) = bb;
     end
    end
    data_beta=data_beta1; i_index=1;
    save('CS2_mNBI_final_run.mat');
  
end

     



