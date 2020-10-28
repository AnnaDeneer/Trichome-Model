function D = diffusionMatrix(ymax,xmax,ynb,xnb,bcond)
% Function for the construction of a diffusion matrix 
%
% Returns a sparse ymax*xmax-by-ymax*xmax coupling matrix for a 
% reaction-diffusion problem in two dimensions.
% xnb and ynb are vectors indicating the connectivity to neighbours in the
% x and y direction, respectively, relative to a reference cell at
% postition (y,x)=(0,0). Bcond indicates the boundary condition, i.e.,
% bcond=0 (periodic boundary condition), bcond=1 (zero-flux boundary
% condition).
%
% It is assumed that the net exchange of a cell is zero, i.e.,
% sum(D) = zeros(1,ymax*xmax) and D(idx,idx) = -1.*sum(D(idx,:)).
%
% Examples:
%
% A reaction-diffusion problem on a two dimensional domain with 50*30
% cells, zero flux boundary conditions and quadratic cells:
% D = couplingMatrix(50,30,[-1 1 0 0],[0 0 -1 1],1);
% 
% Reference: Digiuni, Simona, et al. "A competitive complex formation
% mechanism underlies trichome patterning on Arabidopsis leaves."
% Molecular Systems Biology 4.1 (2008): 217.

if ~(size(ynb)==size(xnb))
  error('Unequal sizes of ynb and xnb')
end

n = ymax*xmax;
D = sparse(n,n);
idx = reshape(1:n,ymax,xmax);
for s=1:length(ynb)
  sidx = circshift(idx,[ynb(s) xnb(s)]);
  ind = sub2ind([n n],idx,sidx);
  if bcond 
    % zeroflux: remove indices of 'outside' cells
    if ynb(s)>0
      ry = 1:ynb(s);
    elseif ynb(s)<0
      ry = ymax:-1:ymax+ynb(s)+1;
    else
      ry = [];
    end
    if xnb(s)>0
      rx = 1:xnb(s);
    elseif xnb(s)<0
      rx = xmax:-1:xmax+xnb(s)+1;
    else
      rx = [];
    end
    ind(ry,:) = [];
    ind(:,rx) = [];
  end
  D(ind) = 1;
end

D = spdiags(-ones(n,1).*sum(D)',0,D);
