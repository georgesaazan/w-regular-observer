% This functions output the possible values of switching signal given its
% last value.
% This function takes an argument i\in [1,...,m] and outputs an array of
% admissible values.
% ex: adm(2).
function a=adm(i)
if i==1
    a=[1,2,3,5];
elseif i==2
    a=[1,2,4,6];
elseif i==3
    a=[1,3,4,7];
elseif i==4
    a=[2,3,4,8];
elseif i==5
    a=[1,5,6,7];
elseif i==6
    a=[2,5,6,8];
elseif i==7
    a=[3,5,7,8];
elseif i==8
    a=[4,6,7,8];
end
end
