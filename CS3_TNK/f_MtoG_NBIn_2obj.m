
function [objs, n_groups, properties, normalNew, tWall, tCPU] ...
        = f_MtoG_NBIn_2obj(weight, beta1, normal1, phi1, fo1, icut1,dim,n_dim, n_sobol,iPoints)
   
     % initial points based on sobol' sequence
     
     data_x_i=[ ...
    0.589049 0.981747;0.981747 0.589049;0.883573 0.883573;0.564505 0.957204;...
    0.957204 0.564505;1.10447 0.220893;0.220893 1.10447;0.773126 0.79767;...
    1.16583	0.404971;0.871301 1.09219;1.09219	0.871301;0.404971 1.16583;...
    0.79767 0.773126;1.1965 0.386563;1.00016 0.975612;0.509282 1.07379;...
    0.901981 0.681087;0.681087 0.901981;1.07379 0.509282;0.386563 1.1965;...
    0.975612 1.00016;1.01243 0.325204];
    data_x_i(end-1,:)=iPoints(1,:);
    data_x_i(end,:)=iPoints(2,:);

    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[]; fail_i=[]; normalNew_i=[];
    
for id=1:n_sobol
    
    data_initial = data_x_i(id,:);
%==========================================================================
% (1) Assign values from Matlab
%==========================================================================
    if nargin == 10
       data_w = weight;
       data_beta = beta1;
       data_norm = normal1;
       data_phi = phi1;
       data_fo = fo1;
       data_icut = icut1;
       data_f_lo = min(iPoints);
       data_f_up = max(iPoints);

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
    %icut    = GAMS.param('yv', data_icut,[c.uels j.uels]);
    x_i     = GAMS.param('x_i',data_initial,j.uels);
    f_u     = GAMS.param('f_u',data_f_up,i.uels);
    f_l     = GAMS.param('f_l',data_f_lo,i.uels);
    
    % write to GDX file
    GAMS.putGDX('inputcs3_nbi_n.gdx',i,j,m,c,beta,normal,f_o,phi,x_i,f_u,f_l);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    %g = GAMS(struct('model','NBI_CS3_WS.gms'));
    g = GAMS(struct('model','NBI_CS3_normal.gms'));
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (3) Write GAMS gdx file that interfaces with MATLAB
%==========================================================================

    % read result variable x if run successful
    if g.status == 0
        %x = GAMS.getGDX('result.gdx','x');
        %x = GAMS.rectify(x, i.uels);

        % Read from text and store
        status1 = copyfile('NBI_CS3_result_n.txt', 'results2.txt');
        %status2 = copyfile('NBI_CS3_WS.log', 'log2.txt');
        status2 = copyfile('NBI_CS3_normal.log', 'log2.txt');
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
         n_data = 18;  %(no.constraint + no.extra)
         data_type = '%f %f %f %f %f %f %f %f %f ';
         store =[];
         for j=1:2
         data_result = textscan(data_resulto, data_type, 9, 'Delimiter', ':');
         store2 = cell2mat(data_result); 
         store(end+1,:)=store2(1,:);
         data_resulto = fgetl(fid);
         end
         
         store = reshape(store',[1,n_data]); 
         store(:,end)=[];
         fclose(fid); 
         fclose(fid2);
         delete('results2.txt');
         delete('log2.txt');
         delete('NBI_CS3_result_n.txt');


         failtest = store(1,end);
         fail_i(end+1,:) = failtest;
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
            continue
         end      
         
         objs_i(end+1,:)       = store(1,1:2);
         n_groups_i(end+1,:)   = store(1,3:4);
         properties_i(end+1,:) = [store(1,5:end-4)];
         normalNew_i(end+1,:) = [store(1,end-3), -1-store(1,end-3)];
         tCPU_i(end+1,:) = store(1,end-2);
         tWall_i(end+1,:) = cputime - tWall_initial;    
               

        end
     

    end
end
        [M, index] = max(properties_i(:,1));
         objs       = objs_i(index,:);
         n_groups   = n_groups_i(index,:);
         properties = properties_i(index,:); 
         normalNew  = normalNew_i(index,:);
         tCPU       = sum(tCPU_i,1);
         tWall      = sum(tWall_i,1);
       
end

