% Function to locate boundary indices for each grid level

function boundary=setboundaryindex_2d(domain,option)

for i=1:option.grids
    

    [I,J]=ndgrid(1:domain(i).N(1),1:domain(i).N(2));
    
    % Find indices of boundary
    boundary(i).index{1}=sub2ind(domain(i).N,I(1,:),J(1,:)); % Top boundary of matrix (x(1))
    boundary(i).index{2}=sub2ind(domain(i).N,I(end,:),J(end,:)); % Bottom boundary of matrix (x(end))
    boundary(i).index{3}=sub2ind(domain(i).N,I(:,1),J(:,1)); % Left boundary of matrix (y(1))
    boundary(i).index{4}=sub2ind(domain(i).N,I(:,end),J(:,end)); % Right boundary of matrix (y(end))
    
end

end