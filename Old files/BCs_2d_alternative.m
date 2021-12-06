% Alternative code to set boundary conditions on 2D matrix in matrix solve
% 
% Sets boundary conditions after kron (not before)
%
% Code may prove useful if ever extended to mixed derivative terms ...



% Dirichlet conditions, can be easily extended to Neumann

% Top
A(1:Nx:(Ny-1)*Nx+1,:)=sparse(Ny,Nx*Ny); % Zero specific rows to BCs
index=sub2ind(size(A),1:Nx:(Ny-1)*Nx+1,1:Nx:(Ny-1)*Nx+1); % Get index
A(index)=1; % Set main diag element of rows to 1

% Bottom
A(Nx:Nx:Nx*Ny,:)=sparse(Ny,Nx*Ny);
index=sub2ind(size(A),Nx:Nx:Nx*Ny,Nx:Nx:Nx*Ny);
A(index)=1;


% Alternatively


[I,J]=ndgrid(1:domain.N(1),1:domain.N(2));
index=cell(length(domain.BCflag),1);

% Find indices of boundary
index{1}=sub2ind(domain.N,I(1,:),J(1,:)); % Top boundary of matrix (x(1))
index{2}=sub2ind(domain.N,I(end,:),J(end,:)); % Bottom boundary of matrix (x(end))
index{3}=sub2ind(domain.N,I(:,1),J(:,1)); % Left boundary of matrix (y(1))
index{4}=sub2ind(domain.N,I(:,end),J(:,end)); % Right boundary of matrix (y(end))

% For spectral matrix ... but not generalised ...
for i=1:length(domain.BCflag)
    
    switch domain.BCflag(i)
        
        case 1 % Dirichlet
            
            A(index{i},:)=0; % Zero specific rows to BCs
            diag_index=sub2ind(size(A),index{i},index{i}); % Get index
            A(diag_index)=1; % Set main diag element of rows to 1
            
        case 2 % Neumann only left :(((
            
            if rem(i,2)~=0 % Is not even?
                A(index{i},:)=0; % Zero specific rows to BCs
                
                [I,J]=ndgrid(index{i},index{3});
                K=I+J-1;
                
                for j=1:Ny
                    
                    A(index{1}(j),K(j,:))=sum(fct(eye(length(index{3}))).*k{i}.^2);
                    
                end

                
            else
%                 A(index{i},:)=sum((-1).^k{1}.*fct(eye(N)).*k{1}.^2); % for x=-1
            end
            
            
        case 3 % Fourier

            % do nothing :)
    end
    
end