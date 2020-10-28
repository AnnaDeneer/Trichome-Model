function J = jMatrix(NVar,D)
% Returns a sparse matrix indicating the pattern of the jacobian matrix.
%
% NVar is the number of species and D is
% the diffusion matrix. The optional input matrix SYSMAT indicates
% the system coupling, i.e., structure of the jacobian matrix of the
% single-cell system. 
%
% Reference: Digiuni, Simona, et al. "A competitive complex formation
% mechanism underlies trichome patterning on Arabidopsis leaves."
% Molecular Systems Biology 4.1 (2008): 217.

sysmat = ones(NVar);
diffmat = speye(NVar);

Js = kron(speye(size(D)),sysmat); % system coupling
Jd = kron(spones(D),diffmat); % diffusive coupling
J = spones(Js+Jd);

  
  
