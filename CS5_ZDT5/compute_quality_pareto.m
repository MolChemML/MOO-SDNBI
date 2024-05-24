function [DM] = compute_quality_pareto(ParetoPoints, f_nadir, f_utopia)
% [crowding_distances] = compute_crowding_distances(F)
%
% Computes crowding distances given a set of l objective function 
% value vectors. The implementation follows the description of:
%
% 'K. Deb, A. Pratap, S. Argawal, and T. Meyarivan. A Fast and
%  Elitist Multi-Objective Genetic Algorithm: NSGA-II. KanGAL
%  Report No. 2000001, 2001.'
%
% IMPORTANT:
% This function assumes that all solutions of ParetoPoints are non-dominated!
%
% Input:
% - ParetoPoints		 - A matrix of P x M, where M is the number
%						  of objectives, andPl is the number of
%						  objective function value vectors of the
%						  solutions.
%
% Output:
% - crowding_distances	- a vector of 1 x l with the crowding distances.  
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011


    [P, M] = size(ParetoPoints);
  
	if (nargin < 3)
          f_nadir =  max(ParetoPoints); % ub
          f_utopia = min(ParetoPoints); % lb
    end
 
    f_ideal_nadir = abs(f_nadir - f_utopia);
    R = max(ParetoPoints) - min(ParetoPoints);
    
    Objs_index = [1:P]';
    Objs_index(:,end+1:end+2) = ParetoPoints;
    
    d_e = zeros(P-1,M);
    for i=1:M
        sort_obj = sortrows(Objs_index,1+i,'ascend');
        for j=1:P-1
            d_e(j,i)= sort_obj(j+1,i+1)-sort_obj(j,i+1);
        end
        mu(i) = sum(d_e(:,i))/(P-1); % mean
        dev(i) = std(d_e(:,i)); % standard deviation
        const_a(i) = dev(i)/mu(i)*f_ideal_nadir(i)/R(i);
    end
      
       DM = 1/P * sum(const_a); % Distribution Measure: The smaller, the better.
    
   



%{
	[M, l] = size(F);
	crowding_distances2 = zeros(1, l);
	for j = 1 : M
		[sort_F, sort_index] = sort(F(j,:));
		%crowding_distances(sort_index(1)) = 1e5;
		%crowding_distances(sort_index(end)) = 1e5;
		for k = 1 : l - 1
			crowding_distances2(sort_index(k)) = crowding_distances2(sort_index(k)) + (F(j, sort_index(k+1)) - F(j, sort_index(k)));
        end
    end
    
    crowding_distances = sum(crowding_distances2,2);
%}
end
