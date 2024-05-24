% Function Description:
% This function interfaces with GAMS to perform optimization tasks. It 
% passes the input data from MATLAB to a GAMS model, runs the GAMS file, 
% and extracts the final values from the results.

% Inputs:
% - weight: A vector representing the weight of objectives.
% - beta1: A vector for the beta parameter.
% - normal1: A vector for the normal parameter.
% - phi1: A matrix for the phi parameter.
% - fo1: A vector for the fo parameter.
% - fok: A matrix for the fok parameter.
% - icut1: A vector for the icut parameter.
% - dim: Dimension of the problem.
% - n_dim: Number of dimensions.
% - n_sobol: Number of Sobol points.

% Outputs:
% - objs: Objectives after optimization.
% - n_groups: Number of groups in the optimized result.
% - properties: Properties of the optimized result.
% - tWall: Total wall clock time.
% - tCPU: Total CPU time.

function [objs, n_groups, properties, tWall, tCPU] ...
        = f_MtoG_mNBIn_2obj(weight, beta1, normal1, phi1, fo1, fok, icut1, dim,n_dim, n_sobol, x_ini, flag_b)

   % 1) initial points only for x1


  n_initial_pt = n_sobol;
 
  
  data_x_i=randi([0,5],n_sobol-1,10);
  data_x_i(end+1,:)= [1 1 1 1 1 1 1 1 1 1 ];
  data_x_i(n_sobol-1,:)= x_ini(1,2:end);
  data_x_i(n_sobol,:)= x_ini(2,2:end);
  
  data_xo_i= randi([0,30],n_sobol,1); 
  data_xo_i(n_sobol-1,:)= x_ini(1,1);
  data_xo_i(n_sobol,:)= x_ini(2,1);
   
    % 2) initial points based on sobol' sequence
    %{
 n_initial_pt = n_sobol;
 fid  = fopen('sobol_results.txt','r');    
 formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
 N=n_initial_pt;
 sobol_initial =  textscan(fid, formatSpec, N,'Delimiter', ':');
 
 data_x_i =cell2mat(sobol_initial);
 fclose(fid);
%}
    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[]; normalNew_i=[];
    
for id=1:n_sobol
    
    data_initial = data_x_i(id,:);
    data_yo_i    = data_xo_i(id,:);
%==========================================================================
% (1) Assign values from input arguments
%==========================================================================
    if nargin >10
       data_w = weight;
       data_beta = beta1;
       data_norm = normal1;
       data_phi = phi1;
       data_fo = fo1;
       data_icut = icut1;
       data_fk = fok;

    else
        disp("Not Enough Input Arguments!")
        break
    end

%==========================================================================
% (2) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================
  tWall_initial = cputime;

    % sets
    i    = GAMS.set('i', {'i1','i2'});
    j    = GAMS.set('j', {'j1','j2','j3','j4','j5','j6','j7','j8','j9','j10'});
    m    = GAMS.set('m', 1:size(data_beta,2));
    %c    = GAMS.set('c', 1:size(data_icut,1));

    % parameters: ex. P = GAMS.param(name, vals, onsets, form)
    beta    = GAMS.param('beta',data_beta,m.uels);
    %w       = GAMS.param('w', data_w,i.uels);
    normal  = GAMS.param('normal', data_norm,i.uels);
    f_o     = GAMS.param('f_o',data_fo,i.uels);
    f_k     = GAMS.param('f_k',data_fk,[i.uels m.uels]);
    phi     = GAMS.param('phi',data_phi, [i.uels m.uels]);
    y_i     = GAMS.param('y_i',data_initial,j.uels);
    yo_i    = GAMS.param('yo_i',data_yo_i); 
    
    % write to GDX file
    GAMS.putGDX('inputcs5_nbi_n.gdx',i,j,m,beta,normal,f_o,f_k, phi,y_i,yo_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    %g = GAMS(struct('model','NBI_CS5_ZDT5_normal.gms'));
    
    if flag_b ==3;
       g = GAMS(struct('model','mNBIn_a_CS5_ZDT5.gms'));
    elseif flag_b ==4;
       g = GAMS(struct('model','mNBIn_b_CS5_ZDT5.gms'));
    end
    
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (4) Read result if GAMS run successful
%==========================================================================

    % read result variable x if run successful
    if g.status == 0

        % Read from text and store
        status1 = copyfile('NBI_CS5_result_n.txt', 'results2.txt');
        if flag_b ==3 ;
            status2 = copyfile('mNBIn_a_CS5_ZDT5.log', 'log2.txt');
        elseif flag_b == 4;
            status2 = copyfile('mNBIn_b_CS5_ZDT5.log', 'log2.txt');
        end 
       
        if status1 == 1
         fid  = fopen('results2.txt','r');
         
         buffer = fgetl(fid);
         data_resulto = fgetl(fid);
         
         n_data = 21;  %(no.constraint + no.extra)
         data_type = '%f %f %f %f %f %f %f ';
         store =[];
         for j=1:3
         data_result = textscan(data_resulto, data_type, 3, 'Delimiter', ':');
         store2 = cell2mat(data_result); 
         store(end+1,:)=store2(1,:);
         data_resulto = fgetl(fid);
         end
         
         store = reshape(store',[1,n_data]); 
         %store(:,end-3:end)=[];
         fclose(fid); 
         delete('results2.txt');
         delete('NBI_CS5_result_n.txt');
         
         failtest = store(1,end);
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
             tCPU_i(end+1,:) = store(1,end-1);
             tWall_i(end+1,:) = cputime - tWall_initial; 
             continue
         end  
         
         if store(1,1)>= min(fok(:,1)) && store(1,1)<= max(fok(:,1)) && store(1,2)>= min(fok(:,2)) && store(1,2)<= max(fok(:,2))
         objs_i(end+1,:)       = store(1,1:dim);
         n_groups_i(end+1,:)   = store(1,dim+1:dim+n_dim);
         properties_i(end+1,:) = [store(1,dim+n_dim+1:end-3)];
         end
         tCPU_i(end+1,:) = store(1,end-1);
         tWall_i(end+1,:) = cputime - tWall_initial;    
         
              
         
        end
      

    end
end
%==========================================================================
% (5) Find maximum properties and corresponding values
%==========================================================================
       if isempty(properties_i)
           objs= zeros([1 dim]);
           n_groups = zeros([1 n_dim]);
           properties = zeros([1 5]);
           tCPU       = sum(tCPU_i,1);
           tWall      = sum(tWall_i,1);
       else
         [M, index] = max(properties_i(:,1));
         objs       = objs_i(index,:);
         n_groups   = n_groups_i(index,:);
         properties = properties_i(index,:); 
         tCPU       = sum(tCPU_i,1);
         tWall      = sum(tWall_i,1);
       end
       
end

