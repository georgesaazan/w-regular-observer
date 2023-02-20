% This function output a matrix that contains minimal observable sequences
% This function takes into its arguments two matrices a and b, the matrix b
% is the matrix whose rows are the observable sequences for which we check
% if they have prefixes in a.
% This function removes the non-minimal observable sequences in b.
% ex: remov([1,2],[1,2,3])
function out=remov(a,b)
n=size(a,2);
out=b; %initialize the matrix for whom we shall remove non-minimal observable sequences
I=[]; %initialize the matrix containing the indexes of the non-minimal observable sequences
for i=1:size(b,1)
   if  ismember(b(i,1:n),a,'rows') %for each sequence in b check if there is a matching prefix in a (non-minimality)
       I=[I,i]; %note the indices of non-minimal sequences in b
   end
end
out(I,:)=[]; %build a new matrix containing the minimal observable sequences in b
end
