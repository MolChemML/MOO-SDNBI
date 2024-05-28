function [index_outer_corrected] = f_outer_test(PPoints,index_outer,A,b)
    
    %pt_test = sortrows(PPoints(index_outer(2:end),3:4),1,'ascend');
    pt_test = PPoints(index_outer(2:end),3:4); 
    
    c_outer = index_outer(1); % if -3 all b's should be >=0  if -9 all b's should be <= 0
    multiplier = 1; index_del=[];
    if c_outer == -9
        multiplier = -1;
    end
    for i=1:(length(index_outer)-1)
         w_correct = A(index_outer(i+1),1:2);
         b_t = b(index_outer(i+1),1);
         w_test = repmat(w_correct,size(pt_test,1),1);
         b_test =  sum(w_test.*pt_test,2) - repmat(b_t,size(pt_test,1),1);
                       
         if c_outer == -3
             if prod(b_test >= 0) == 0
                 index_del(end+1)= i+1;
             end
         elseif c_outer == -9
             if prod(b_test <= 0) == 0
                 index_del(end+1)= i+1;
             end
         end
    end
    index_outer(index_del)=[];
    index_outer_corrected = index_outer;
end