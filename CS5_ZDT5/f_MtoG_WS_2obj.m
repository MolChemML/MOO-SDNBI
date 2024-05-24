% gams-matlab interface:  "NBI_3obj_int.gms" (for furter info.:
% https://www.youtube.com/watch?v=s755kr8MH0Y)
% A gdx file only accept structure format.
% system 'gams NBI_3obj_int lo=3 gdx=tmp'
% gdxInfo tmp

function [objs, n_groups, properties, tWall, tCPU] ...
        = f_MtoG_WS_2obj(weight,dim,n_dim, n_sobol)
 
  n_initial_pt = n_sobol;
 
  data_x_i=randi([0,5],n_sobol-2,10);
  data_x_i(end+1,:)= [1 1 1 1 1 1 1 1 1 1 ];
  data_x_i(end+1,:)= [1 1 1 1 1 1 1 1 1 1 ];
  %data_x_i= repmat([1],n_sobol,10);
  data_xo_i= randi([0,30],n_sobol,1);

  
 
  %{
 fid  = fopen('sobol_results.txt','r');    
 formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
 N=n_initial_pt;
 sobol_initial =  textscan(fid, formatSpec, N,'Delimiter', ':');
 
 data_x_i =cell2mat(sobol_initial);
 fclose(fid);
 %}
 
 objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[]; fail_i=[];
%==========================================================================
% (1) Assign values from Matlab
%==========================================================================
for i=1:n_initial_pt
    
    data_initial = data_x_i(i,:);
    data_yo_i = data_xo_i(i,:);
    
    if nargin > 3
       data_w = weight;

    else 
        data_w = [1 0 ];    % default value
    end
    
    data_beta = [0.333 0.333]; % default value  
    data_norm = [-0.7071, -0.7071]; % default value
    data_fo = [-6, -6];
    data_phi = [ 6  0  ;   0  1];

%==========================================================================
% (2) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================
  tWall_initial = cputime;
    % sets
    i    = GAMS.set('i', {'i1','i2'});
    j    = GAMS.set('j', {'j1','j2','j3','j4','j5','j6','j7','j8','j9','j10'});
    k    = GAMS.set('k', {'k1','k2'});

    % parameters: ex. P = GAMS.param(name, vals, onsets, form)
    beta    = GAMS.param('beta',data_beta,i.uels);
    w       = GAMS.param('w', data_w,i.uels);
    normal  = GAMS.param('normal', data_norm,i.uels);
    f_o     = GAMS.param('f_o',data_fo,i.uels);
    phi     = GAMS.param('phi',data_phi, [i.uels j.uels]);
    y_i     = GAMS.param('y_i',data_initial,j.uels);
    yo_i    = GAMS.param('yo_i',data_yo_i);

    % write to GDX file
    GAMS.putGDX('inputcs5_ws.gdx',i,j,beta,w,normal,f_o,phi,y_i,yo_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    g = GAMS(struct('model','NBI_CS5_ZDT5_WS.gms'));
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (3) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================

    % read result variable x if run successful
    if g.status == 0
        %x = GAMS.getGDX('result.gdx','x');
        %x = GAMS.rectify(x, i.uels);

        % Read from text and store
        status1 = copyfile('NBI_CS5_WS_result.txt', 'results2.txt');
        status2 = copyfile('NBI_CS5_ZDT5_WS.log', 'log2.txt');
        if status1 == 1
         fid  = fopen('results2.txt','r');
         fid2 = fopen('log2.txt','r');
         
         buffer = fgetl(fid);
         data_resulto = fgetl(fid);
         cpu = textscan(fid2, '%q %q', 'Delimiter', ':');
         
          %'Wall clock time' = 
          %{
          str_compare = ['Total CPU time used'];
          for j=1:size(cpu{1,1}(:,1),1)
               if  strncmp(str_compare, cell2mat(cpu{1,1}(j,1)), 12);
                   tCPU = str2double(cell2mat(cellstr(cpu{1,2}(j,1))));
                   break
               end
          end
         tCPU=0.1;
          %}
         n_data =21;  %(no.constraint + no.extra)
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
         delete('NBI_CS5_WS_result.txt');
         
         
          failtest = store(1,end);
          fail_i(end+1,:)=failtest;
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
       [M, index] = min(properties_i(:,1));
         objs       = objs_i(index,:);
         n_groups   = n_groups_i(index,:);
         properties = properties_i(index,:); 
         tCPU       = sum(tCPU_i,1);
         tWall      = sum(tWall_i,1);
       end
end

%{
%==========================================================================
% sets
t      = GAMS.set('t', 1:8760);
i      = GAMS.set('i', {'pv', 'windon', 'windoff'});

% parameters
cs     = GAMS.param('cs',100); % cost of storage (EUR/MWh)
c      = GAMS.param('c',[3000 1500 2500],i.uels); % cost of plant (EUR/MWh)
demand = GAMS.param('demand',rand(8760,1),t.uels);

% renewable timeseries
values = [ ... 
    min(max(0, sin((1:8760)'/24*3.14/2).^4+0.15*randn(8760,1)), 1), ...
    min(max(0, rand(8760,1)), 1), ...
    min(max(0, rand(8760,1).^0.25), 1) ];
onset = [ t.uels i.uels ];
cf = GAMS.param('cf', values, onset);
clear values onset;

% write to GDX file
GAMS.putGDX('input.gdx',t,i,c,cs,demand,cf);

% run GAMS model
g = GAMS(struct('model','fuelstation.gms'));
g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"

% read result variable x if run successful
if g.status == 0
	x = GAMS.getGDX('result.gdx','x');
	x = GAMS.rectify(x, i.uels);

    % create bar chart of installed plant capacities
	bar(1000*x.val);
	set(gca,'XTickLabel',x.uels{1});
	ylabel('Installed capacity (kW)');
end
%}