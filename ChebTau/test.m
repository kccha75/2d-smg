BB=zeros(size(AA));
BB(1,:)=1;
BB=sparse(BB);

CC=AA+BB;