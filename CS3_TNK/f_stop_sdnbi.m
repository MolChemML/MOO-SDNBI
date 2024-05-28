function flag_stop = ...
       f_stop_sdnbi(facetIndices, dMax, sMax, dMeasure, dTol, sTol, data_subPrbs, iter, n)


      
      % for the subregion that satisfy the tolerance criteria (inner and
      % outer distance, dMax, and distribution, sMAX), opt out from the
      % consideration
      flag_stop = 0;
      if (dMax < dMeasure || sMax < sTol)
         if isempty(facetIndices)
            sMax =0;
         end
         if (sMax < sTol)
            %dMax= dTol+1; 
            idx_test=sum(~cellfun(@isempty,data_subPrbs),2); 
            if  size(idx_test,1) >= iter
                flag_stop = 1;
                return;
            elseif (n == idx_test(iter-1) && size(idx_test,1) < iter)
                    %dMax = 0;
                    flag_stop = 2;
                    return;
            else       
                %dMax= dTol+1;
                flag_stop = 1;
                return;
            end
         end
         %dMax= dTol+1;
         flag_stop=3;
      end
      

      
end
