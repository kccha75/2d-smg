% function to locate maximum and local wiggle minimum for DJL tabletop solution
function [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v,minindex_x,maxindex_x)

[maxi,index_z]=max(max(v)); % locate maximum, and index for max (for z)

if isempty(maxindex_x) && isempty(minindex_x)

    % locate local minimum at slice v1(:,indexmax)
    % looks at middle +-50 for minimum
    [mini,minindex_x]=min(v(size(v,1)/2-50:size(v,1)/2+50,index_z)); 
    % adjust location
    minindex_x=size(v,1)/2-50+minindex_x-1;

    % locate local minimum at slice v1(:,indexmax)]
    [maxi,maxindex_x]=max(v(size(v,1)/2-50:size(v,1)/2+50,index_z)); 
    % adjust location
    maxindex_x=size(v,1)/2-50+maxindex_x-1;

else

    [mini]=v(minindex_x,index_z); % locate local minimum at slice v1(:,indexmax)
    [maxi]=v(maxindex_x,index_z); % locate local minimum at slice v1(:,indexmax)]

end

end