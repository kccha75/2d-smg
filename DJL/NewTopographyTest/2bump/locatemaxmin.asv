% function to locate maximum and local wiggle minimum for DJL tabletop solution
%
% Note that local minimum and maximum location input is to make sure secant
% method can return a negative number for difference=maxi-mini
%
% This assumes the local minimum and maximum location stays fixed, and that
% for some value of delta can switch minimum to maximum and vice versa at
% those locations (this seems to be true)
%
% Inputs:
% v - DJL tabletop solution
% minindex_x - x index for local minimum (optional input)
% maxindex_x - x index for local maximum (optional input)
%
% Outputs:
% mini - local minimum value
% maxi - local maximum value
% minindex_x - x index for local minimum (if inputted, output is same)
% maxindex_x - x index for local maximum (if inputted, output is same)
% index_z - z index where local minimum and maximum take place

function [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v,minindex_x,maxindex_x)

[maxi,index_z]=max(max(v)); % locate maximum, and index for max (for z)

% if empty inputs ... find the indices!
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

    [mini]=v(minindex_x,index_z); % locate local minimum at slice v(:,indexmax)
    [maxi]=v(maxindex_x,index_z); % locate local minimum at slice v(:,indexmax)

end

end