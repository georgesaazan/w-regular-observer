function out =mb(s,A,C,rhoo)
s=flip(s);
out=[];
for i=1:length(s)
    out=[out;mseqb(s,A,C,rhoo)];
    s(1)=[];
end
out=flip(out)'*flip(out);
end
function out=mseqb(s,A,C,rhoo)
out=C{s(1)};
if length(s)~=1
for i=1:length(s)-1
    out=out*A{s(i+1)}/rhoo;
end
end
end
