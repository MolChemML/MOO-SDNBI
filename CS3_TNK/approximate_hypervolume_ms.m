function [hv] = approximate_hypervolume_ms(F, ub, samples)
% [hv] = approximate_hypervolume_ms(F, ub, samples)
%
% Computes the hypervolume (or Lebesgue measure) of the M x P
% matrix F of P vectors of M objective function values by means
% of a Monte-Carlo approximation method.
%
% IMPORTANT:
%   Considers Minimization of the objective function values!
%
% Input:
% - F              - An M x P matrix where each of the P columns 
%                    represents a vector of M objective function values
%                    (note that this function assumes that all 
%                    solutions of F are non-dominated).
% - ub             - Optional: Upper bound reference point (default:
%                    the bouFndary point containing the maximum of F
%                    for each objective).
% - samples        - Optional: The number of samples used for the Monte-
%                    Carlo approximation (default: 100000).
%
% Output:
% - hv             - The hypervolume (or lebesgue measure) of F.
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	if (nargin < 2), ub = (max(F,[],2)); end
	if (nargin < 3), samples = 10000; end
    %ub = [0.9 1 1 1 1]';   % fixed referecne point: LL !!!!!!!
	[M, P] = size(F);
	samples = 100000;
	lb = min(F')';

    if M == 2 % for biobjective function, easy to calucate hypervolume.
       objs = F';
       objs_sort = sortrows(objs,1,'ascend');
       hv_exact = 0;
       counter = 1;
       ref_pt = ub;
       
       hv_exact = (ref_pt(1)-objs_sort(1,1))*(ref_pt(2)-objs_sort(1,2));
       for i=1:P-1
           fragment = (ref_pt(1)-objs_sort(i+1,1))*(objs_sort(i,2)-objs_sort(i+1,2));
           hv_exact = hv_exact + fragment;
           counter = counter + 1;
       end
       hv = hv_exact;
    else 
       F_samples = repmat(lb,1,samples) + rand(M,samples) .* repmat((ub - lb),1,samples);
	   is_dominated_count = 0;
	   for i = 1:samples
            for j = 1:P
                if (dominates(F(:,j), F_samples(:,i)))
                    is_dominated_count = is_dominated_count + 1;
                    break;
                end
            end
       end

	   hv = prod(ub - lb) * (is_dominated_count / samples);
    end

end
