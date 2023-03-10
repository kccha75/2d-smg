% function to locate maximum and local wiggle minimum for DJL tabletop solution
function difference=locatemaxmin(v)

[maxi,indexmax]=max(max(v)); % locate maximum, and index for max (for z)

mini=min(v(size(v,1)/2-50:size(v,1)/2+50,indexmax)); % locate local minimum at slice v1(:,indexmax)

difference=maxi-mini;

end