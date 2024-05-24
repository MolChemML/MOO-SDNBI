% gams-matlab interface:  "NBI_3obj_int.gms" (for furter info.:
% https://www.youtube.com/watch?v=s755kr8MH0Y)
% A gdx file only accept structure format.
% system 'gams NBI_3obj_int lo=3 gdx=tmp'
% gdxInfo tmp

function [objs, n_groups, properties, tWall, tCPU] ...
        = f_MtoG_WS_2obj(weight, n_sobol)

    data_x_i=[2.5,1.5;1.25,2.25;1.875,1.125;3.125,0.375;0.625,1.875 ; ...
        0.9375,0.9375;2.1875,1.6875;1.5625,0.5625;2.8125,1.3125;0.3125,2.8125; ...
        0.46875,1.40625;1.71875,2.15625;2.34375,0.28125;1.09375,2.53125;0.78125,0.46875; ...
         1.40625,0.84375;2.65625,0.09375;0.15625,1.59375;0.234375,0.796875;...
1.48438,1.54688;2.10938,0.421875;3.35938,1.17188;0.859375,2.67188;...
1.17188,0.234375;1.79688,1.35938;3.04688,0.609375;0.546875,2.10938;...
0.390625,0.703125;2.26563,1.07813;3.51563,0.328125;1.01563,1.82813;0.703125,1.26563;...
1.95313,2.01563	;1.32813,0.140625;2.57813,0.890625;0.078125,2.39063;0.117188,1.19531;0,0];	


    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[]; fail_i=[];
    
 for i=1:n_sobol
     data_initial = data_x_i(i,:);
%==========================================================================
% (1) Assign values from Matlab
%==========================================================================
    if nargin > 1
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
    j    = GAMS.set('j', {'j1','j2'});
    k    = GAMS.set('k', {'k1','k2'});

    % parameters: ex. P = GAMS.param(name, vals, onsets, form)
    beta    = GAMS.param('beta',data_beta,i.uels);
    w       = GAMS.param('w', data_w,i.uels);
    normal  = GAMS.param('normal', data_norm,i.uels);
    f_o     = GAMS.param('f_o',data_fo,i.uels);
    phi     = GAMS.param('phi',data_phi, [i.uels j.uels]);
    x_i     = GAMS.param('x_i',data_initial,j.uels);
    
    % write to GDX file
    GAMS.putGDX('inputcs2_ws.gdx',i,j,beta,w,normal,f_o,phi,x_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    g = GAMS(struct('model','NBI_CS2_WS.gms'));
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (3) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================

    % read result variable x if run successful
    if g.status == 0
        %x = GAMS.getGDX('result.gdx','x');
        %x = GAMS.rectify(x, i.uels);

        % Read from text and store
        status1 = copyfile('NBI_CS2_WS_result.txt', 'results2.txt');
        status2 = copyfile('NBI_CS2_WS.log', 'log2.txt');
        if status1 == 1
         fid  = fopen('results2.txt','r');
         fid2 = fopen('log2.txt','r');
         
         buffer = fgetl(fid);
         data_resulto = fgetl(fid);
         cpu = textscan(fid2, '%q %q', 'Delimiter', ':');
         
          %'Wall clock time' = 
          str_compare = ['Total CPU time used'];
          for j=1:size(cpu{1,1}(:,1),1)
               if  strncmp(str_compare, cell2mat(cpu{1,1}(j,1)), 12);
                   tCPU = str2double(cell2mat(cellstr(cpu{1,2}(j,1))));
                   break
               end
          end
         tCPU=0.1;
         data_result = textscan(data_resulto, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter', ':');
         store = cell2mat(data_result); 
         fclose(fid); 
         fclose(fid2);
         delete('results2.txt');
         delete('log2.txt');
         delete('NBI_CS2_WS_result.txt');
 
         
           
        
         failtest = store(1,end);
         fail_i(end+1,:) = failtest; 
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
            continue
         end   
                  
         objs_i(end+1,:)       = store(1,1:2);
         n_groups_i(end+1,:)   = store(1,3:4);
         properties_i(end+1,:) = [store(1,5:end-3),0,0]; 
         tCPU_i(end+1,:) = store(1,end-2);
         tWall_i(end+1,:) = cputime - tWall_initial;     
   
       
         

        end
       % create bar chart of installed plant capacities
    end
    
 end
    
         [M, index] = min(properties_i(:,1));
         objs       = objs_i(index,:);
         n_groups   = n_groups_i(index,:);
         properties = properties_i(index,:); 
         tCPU       = sum(tCPU_i);
         tWall      = sum(tWall_i);  
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