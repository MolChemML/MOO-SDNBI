% gams-matlab interface:  "NBI_3obj_int.gms" (for furter info.:
% https://www.youtube.com/watch?v=s755kr8MH0Y)
% A gdx file only accept structure format.
% system 'gams NBI_3obj_int lo=3 gdx=tmp'
% gdxInfo tmp

function [objs, n_groups, properties, tWall, tCPU] ...
        = f_MtoG_oNBI_2obj(weight, beta1, normal1, phi1, fo1, fkl, icut1,dim,n_dim,data_x_i)
    

    %data_x_i=[0.1 0.9; 0.2 0.8; 0.3 0.7; 0.5 0.5; 0.6 0.9; 0.8 0.9; 0.9 0.6; 0.99 0.2];
    %data_x_i=[0.1 0.9; 0.2 0.8; 0.3 0.7; 0.5 0.5; 0.7 0.3; 0.8 0.2; 0.9 0.1];
    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[];
    
for i=1:size(data_x_i,1)
      data_initial = data_x_i(i,:);
%==========================================================================
% (1) Assign values from Matlab
%==========================================================================
    if nargin > 4
       data_w = weight;
       data_beta = beta1;
       data_norm = normal1;
       data_phi = phi1;
       data_fo = fo1;
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
  tWall_initial = cputime;

    % sets
    i    = GAMS.set('i', {'i1','i2'});
    j    = GAMS.set('j', {'j1','j2'});
    m    = GAMS.set('m', 1:size(data_beta,2));
    c    = GAMS.set('c', 1:size(data_icut,1));

    % parameters: ex. P = GAMS.param(name, vals, onsets, form)
    beta    = GAMS.param('beta',data_beta,m.uels);
    %w       = GAMS.param('w', data_w,i.uels);
    normal  = GAMS.param('normal', data_norm,i.uels);
    f_o     = GAMS.param('f_o',data_fo,i.uels);
    phi     = GAMS.param('phi',data_phi, [i.uels m.uels]);
    icut    = GAMS.param('yv', data_icut,[c.uels j.uels]);
    x_i     = GAMS.param('x_i',data_initial,j.uels);

    % write to GDX file
    GAMS.putGDX('inputcs3_onbi.gdx',i,j,m,c,beta,normal,f_o,phi,icut,x_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    %g = GAMS(struct('model','NBI_CS3_WS.gms'));
    g = GAMS(struct('model','NBI_CS3_o.gms'));
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (3) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================

    % read result variable x if run successful
    if g.status == 0
        %x = GAMS.getGDX('result.gdx','x');
        %x = GAMS.rectify(x, i.uels);

        % Read from text and store
        status1 = copyfile('NBI_CS3_o_result.txt', 'results2.txt');
        %status2 = copyfile('NBI_CS3_WS.log', 'log2.txt');
        status2 = copyfile('NBI_CS3_o.log', 'log2.txt');
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
         %tCPU= 0.1; %temporary
         %}
         data_result = textscan(data_resulto, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f  %f', 'Delimiter', ':');
         store = cell2mat(data_result(1,1:end)); 
         fclose(fid); 
         fclose(fid2);
         delete('results2.txt');
         delete('log2.txt');
         delete('NBI_CS3_o_result.txt');
         
         
         failtest = store(1,end);
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
            continue
         end      
         
         objs_i(end+1,:)       = store(1,1:2);
         n_groups_i(end+1,:)   = store(1,3:4);
         properties_i(end+1,:) = [store(1,5:end-3)]; 
         tCPU_i(end+1,:) = store(1,end-2);
         tWall_i(end+1,:) = cputime - tWall_initial;       
            
         

        end
       % create bar chart of installed plant capacities

    end
end
       [M, index] = max(properties_i(:,1));
         objs       = objs_i(index,:);
         n_groups   = n_groups_i(index,:);
         properties = properties_i(index,:); 
         tCPU       = tCPU_i(index,:);
         tWall      = tWall_i(index,:);
       
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