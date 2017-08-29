function [permutation_matrix, sorted_array] = sort_perm_matrix( array, tol )

[sorted_array,permutation_matrix] = sortrows(array,1);

cnt = 1;
old_x = sorted_array(1,1);

for i=2:size(sorted_array,1)
    
    i
     
    if abs(old_x - sorted_array(i,1)) < max( abs(sorted_array(i,1)), old_x)*tol*eps()+tol*realmin()
       
        cnt = cnt+1;
    else % different than previous
        
        %fix previously tied elements
        if ( cnt > 1 ) 
          [~,perm] = sort( sorted_array(i-cnt: i -1, 2) );
          sorted_array(i-cnt:i -1,:) = sorted_array(i-cnt -1 + perm,:);
          permutation_matrix(i-cnt:i-1) = i-cnt -1 + perm;
        end
        cnt = 1;
        old_x = sorted_array(i,1);
    end
    
end

if ( cnt > 1 ) 
   n_entries = length(permutation_matrix);
   [~,perm] = sort( sorted_array(n_entries-cnt+1: n_entries, 2) );
   sorted_array(n_entries-cnt+1:n_entries,:) = sorted_array(n_entries - cnt + perm,:);
   permutation_matrix(n_entries-cnt+1:n_entries) = n_entries - cnt + perm;
end

end