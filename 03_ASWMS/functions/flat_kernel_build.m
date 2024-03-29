function F = flat_kernel_build(xnode, ynode, N)
% Build flattening kernel for tomography problem
% (i.e., smooth first derivative)
% 
% N is the total number of model parameters in the final G matrix
% JBR 7/22/18
%
% xnode : latitude nodes (yes, it appears backwards...)
% ynode : longitude nodes (yes, it appears backwards...)

Nx = length(xnode);
Ny = length(ynode);

% Smoothing in y
[i,j] = ndgrid(1:Nx,1:(Ny-1)); % set up grid of points that correspond to grid size, smaller by 1 than y
ind = j(:) + Ny*(i(:)-1); % the indices of each point
dy = diff(ynode); % grid spacing
% dy1 = dy(j(:)-1);
dy2 = dy(j(:)); % grid spacing mapped to each adjacent pair of points
Areg = sparse(repmat(ind,1,2),[ind,ind+1], ...
    [-1./dy2, 1./dy2],N,N);
% Smoothing in x
[i,j] = ndgrid(1:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
dx = diff(xnode);
% dx1 = dx(i(:)-1);
dx2 = dx(i(:));
Areg = [Areg;sparse(repmat(ind,1,2),[ind,ind+Ny], ...
    [-1./dx2, 1./dx2],N,N)];
% F=Areg;

F=sparse(Nx*Ny*2*2,Nx*Ny*2);
for n=1:size(Areg,1)
    ind=find(Areg(n,:)~=0);
    F(2*n-1,2*ind-1)=Areg(n,ind);
    F(2*n,2*ind)=Areg(n,ind);
end

if 0
    for ieq = 1:size(F,1)
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                Fgrid(i,j)= F(ieq,n);
            end
        end
        figure(99); clf;
        imagesc(xnode,ynode,Fgrid);
        colorbar
        pause;
    end
end

return
end
