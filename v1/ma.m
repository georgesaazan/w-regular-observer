% This function outputs a matrix
% [C(i_1);C(i_2)*A(i_1);...;C(i_n)*A(i_(n-1))...*A_1] corresponding to a
% sequence (i_1,...,i_n).
% This function takes the following arguments respectively: the vector (i_1,...,i_n), the set of matrices A, the set of output matrices C.
% ex: ma([2,1],{eye(2),ones(2)},{[1,1],[1,0]})
function out =ma(s,A,C)
s=flip(s); %s=(i_1,...,i_n)-->(i_n,...,i_1)
out=[]; %initialize the output matrix
for i=1:length(s)
    out=[out;mseq(s,A,C)];%compute C(i_n)*A(i_(n-1))...*A_1
    s(1)=[]; %remove the first indice in the sequence to obtain (i_(n-1),...,i_1)
end %repeat the process to obtain [C(i_n)*A(i_(n-1))...*A_1;...;C(i_1)]
out=flip(out);% obtain [C(i_1);C(i_2)*A(i_1);...;C(i_n)*A(i_(n-1))...*A_1]
end
function out=mseq(s,A,C)%for s=(i_n,...,i_1) compute C(i_n)*A(i_(n-1))...*A_1
out=C{s(1)};
if length(s)>1
for i=1:length(s)-1
    out=out*A{s(i+1)};%construct C(i_n)*A(i_(n-1))...*A_1 where n=length(s)
end
end
end
