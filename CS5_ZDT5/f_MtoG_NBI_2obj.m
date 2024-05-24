
function [objs, n_groups, properties, tWall, tCPU] ...
        = f_MtoG_NBI_2obj(weight, beta1, normal1, phi1, fo1, fk1, icut1,dim,n_dim, n_sobol)

   % 1) initial points only for x1
  
  n_initial_pt = n_sobol;
 
  data_x_i=randi([0,5],n_sobol-2,10);
  data_x_i(end+1,:)= [1 1 1 1 1 1 1 1 1 1 ];
  data_x_i(end+1,:)= [1 1 1 1 1 1 1 1 1 1 ];
  %data_x_i= repmat([1],n_sobol,10);
  data_xo_i= randi([0,30],n_sobol,1);
  
   %{
     % 2) initial points based on sobol' sequence
 n_initial_pt = n_sobol;
 fid  = fopen('sobol_results.txt','r');    
 formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
 N=n_initial_pt;
 sobol_initial =  textscan(fid, formatSpec, N,'Delimiter', ':');
 
 data_x_i =cell2mat(sobol_initial);
 fclose(fid);
%}  
    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[];
    
for i=1:n_initial_pt
    
    data_initial = data_x_i(i,:);
    data_yo_i = data_xo_i(i,:);
%==========================================================================
% (1) Assign values from input arguments
%==========================================================================
    if nargin > 9
       data_w = weight;
       data_beta = beta1;
       data_norm = normal1;
       data_phi = phi1;
       data_fo = fo1;
       data_fk = fk1;
       data_icut = icut1;

    elseif nargin == 3
        data_w = [1 0]; % default value
        data_beta = beta1;
        data_norm = normal1;
    else 
        data_w = [1 0 ];    % default value
        data_beta = [0.5 0.5]; % default value  
        data_norm = [-0.7071, 0]; % default value
        data_fo = [-6, 0 ];
        data_phi = [ 6  0;  0  1];
    end

%==========================================================================
% (2) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================
   % Measure  time
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
    %icut    = GAMS.param('yv', data_icut,[c.uels j.uels]);
    y_i     = GAMS.param('y_i',data_initial,j.uels);
    yo_i    = GAMS.param('yo_i',data_yo_i);    

    % write to GDX file
    GAMS.putGDX('inputcs5_nbi.gdx',i,j,m,beta,normal,f_o,f_k,phi,y_i,yo_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    g = GAMS(struct('model','NBI_CS5_ZDT5.gms'));
    g.run; 


%==========================================================================
% (4) Read result if GAMS run successful
%==========================================================================

    % read result variable x if run successful
    if g.status == 0
        %x = GAMS.getGDX('result.gdx','x');
        %x = GAMS.rectify(x, i.uels);

        % Read from text and store
        status1 = copyfile('NBI_CS5_result.txt', 'results2.txt');
        status2 = copyfile('NBI_CS5_ZDT5.log', 'log2.txt');
        if status1 == 1
         fid  = fopen('results2.txt','r');
         fid2 = fopen('log2.txt','r');
         
         buffer = fgetl(fid);
         data_resulto = fgetl(fid);
         cpu = textscan(fid2, '%q %q', 'Delimiter', ':');
         

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
         fclose(fid); 
         fclose(fid2);
         delete('results2.txt');
         delete('log2.txt');
         delete('NBI_CS5_result.txt');
         
         failtest = store(1,end);
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
             tCPU_i(end+1,:) = store(1,end-2);
             tWall_i(end+1,:) = cputime - tWall_initial; 
            continue
         end      
         
         objs_i(end+1,:)       = store(1,1:dim);
         n_groups_i(end+1,:)   = store(1,dim+1:dim+n_dim);
         properties_i(end+1,:) = [store(1,dim+n_dim+1:end-3)]; 
         tCPU_i(end+1,:) = store(1,end-2);
         tWall_i(end+1,:) = cputime - tWall_initial;       
         

        end
       % create bar chart of installed plant capacities

    end
end
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
