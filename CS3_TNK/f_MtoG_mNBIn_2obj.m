% gams-matlab interface:  "NBI_3obj_int.gms" (for furter info.: https://www.youtube.com/watch?v=s755kr8MH0Y)
% A gdx file only accept structure format. system 'gams NBI_3obj_int lo=3 gdx=tmp'
% gdxInfo tmp

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
        = f_MtoG_mNBIN_2obj(weight, beta1, normal1, phi1, fo1, fok, icut1, dim, n_dim, n_sobol, flag_b)

    % Define initial data for multi-start
    data_x_i=[ ...
    0.589049 0.981747;0.981747 0.589049;0.883573 0.883573;0.564505 0.957204;...
    0.957204 0.564505;1.10447 0.220893;0.220893 1.10447;0.773126 0.79767;...
    1.16583	0.404971;0.871301 1.09219;1.09219	0.871301;0.404971 1.16583;...
    0.79767 0.773126;1.1965 0.386563;1.00016 0.975612;0.509282 1.07379;...
    0.901981 0.681087;0.681087 0.901981;1.07379 0.509282;0.386563 1.1965;...
    0.975612 1.00016;1.01243 0.325204; ...
    0.325204,1.01243;0.76699,0.76699;0.693359,1.08606;1.08606,0.693359;0.791534,0.791534;1.1873,0.395767; ...
    0.745515,1.0339;1.13821,0.641204;0.843689,0.935728;1.06458,0.273049;0.377359,1.15662;0.819146,0.518485; ...
    0.635068,0.849825;1.02777,0.457126;0.193282,0.997087;1.05231,0.825281;0.561437,1.11981;0.954136,0.727107; ...
    0.217825,0.972544;0.997087,0.193282];

    data_x_i=[0.45,0.93;0,0;2.5,1.5;1.25,2.25;1.875,1.125;3.125,0.375;0.625,1.875 ; ...
        0.9375,0.9375;2.1875,1.6875;1.5625,0.5625;2.8125,1.3125;0.3125,2.8125; ...
        0.46875,1.40625;1.71875,2.15625;2.34375,0.28125;1.09375,2.53125;0.78125,0.46875; ...
         1.40625,0.84375;2.65625,0.09375;0.15625,1.59375;0.234375,0.796875;...
1.48438,1.54688;2.10938,0.421875;3.35938,1.17188;0.859375,2.67188;...
1.17188,0.234375;1.79688,1.35938;3.04688,0.609375;0.546875,2.10938;...
0.390625,0.703125;2.26563,1.07813;3.51563,0.328125;1.01563,1.82813;0.703125,1.26563;...
1.95313,2.01563	;1.32813,0.140625;2.57813,0.890625;0.078125,2.39063;0.117188,1.19531];

    % Initialize output  variables
    objs_i=[]; n_groups_i=[]; properties_i=[]; tWall_i=[]; tCPU_i=[];

    % Iterate through the data 
for i=1:n_sobol
      data_initial = data_x_i(i,:);
%==========================================================================
% (1) Assign values from input arguments
%==========================================================================
    if nargin > 4
       data_w = weight;
       data_beta = beta1;
       data_norm = normal1;
       data_phi = phi1;
       data_fo = fo1;
       data_fk = fok;
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
    % Measure wall clock time
    tWall_initial = cputime;

    % Define GAMS sets
    i    = GAMS.set('i', {'i1','i2'}); 
    j    = GAMS.set('j', {'j1','j2'});
    %m    = GAMS.set('m', 1:size(data_beta,2));
    m    = GAMS.set('m', {'m1','m2'});
    

    % Define GAMS parameters
    beta    = GAMS.param('beta',data_beta,m.uels);
    %w       = GAMS.param('w', data_w,i.uels);
    normal  = GAMS.param('normal', data_norm,i.uels);
    f_o     = GAMS.param('f_o',data_fo,i.uels);
    f_k     = GAMS.param('f_k',data_fk,[i.uels m.uels]);
    phi     = GAMS.param('phi',data_phi, [i.uels m.uels]);
    %icut    = GAMS.param('yv', data_icut,[c.uels j.uels]);
    x_i     = GAMS.param('x_i',data_initial,j.uels);

    % write above to GDX file
    GAMS.putGDX('inputcs3_nbi_n.gdx',i,j,m,beta,normal,f_o,f_k,phi,x_i);


%==========================================================================
% (3) Run GAMS: --.gms
%==========================================================================

    % run GAMS model
    if flag_b == 3;
       g = GAMS(struct('model','NBI_CS3_mNBIn_L.gms'));   
    elseif flag_b == 4;
       g = GAMS(struct('model','NBI_CS3_mNBIn_R.gms'));
    end
    g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"


%==========================================================================
% (4) Read result if GAMS run successful
%==========================================================================


    if g.status == 0
        % Read from text and store
        status1 = copyfile('mNBIn_CS3_result.txt', 'results2.txt');
        %status2 = copyfile('NBI_CS3_WS.log', 'log2.txt');
        if flag_b ==3 ;
            status2 = copyfile('NBI_CS3_mNBIn_L.log', 'log2.txt');
        elseif flag_b == 4;
            status2 = copyfile('NBI_CS3_mNBIn_R.log', 'log2.txt');
        end
        
        if status1 == 1
         fid  = fopen('results2.txt','r');
         fid2 = fopen('log2.txt','r');
         
         buffer = fgetl(fid);
         data_resulto = fgetl(fid);
         cpu = textscan(fid2, '%q %q', 'Delimiter', ':');
         
         data_result = textscan(data_resulto, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ':');
         store = cell2mat(data_result(1,1:end)); 
         fclose(fid); 
         fclose(fid2);
         delete('results2.txt');
         delete('log2.txt');
         delete('mNBIn_CS3_result.txt');
         
         % Check for valid results
         failtest = store(1,end);
         if ~(failtest==1 || failtest==2 || failtest==8 || failtest==15 || failtest==16)
            continue;
         end      
         
         % Store results
         objs_i(end+1,:)       = store(1,1:2);
         n_groups_i(end+1,:)   = store(1,3:4);
         properties_i(end+1,:) = [store(1,5:end-3)]; 
         tCPU_i(end+1,:) = store(1,end-2);
         tWall_i(end+1,:) = cputime - tWall_initial;       
        end
    end
end

%==========================================================================
% (5) Find maximum properties and corresponding values
%==========================================================================
       [M, index] = max(properties_i(:,1));
        objs       = objs_i(index,:);
        n_groups   = n_groups_i(index,:);
        properties = properties_i(index,:) 
        tCPU       = sum(tCPU_i,1);
        tWall      = sum(tWall_i,1);
        

end
