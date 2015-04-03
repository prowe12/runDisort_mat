
function [DTAUC_new, SSALB_new, NPmom_new, Pmom_new] = combine2(...
  DTAUC_gas, SSALB_gas, NPmom_gas, Pmom_gas, ...
  DTAUC_wat, SSALB_wat, NPmom_wat, Pmom_wat)
%
% Inputs:
%   DTAUC_gas      vector, Nlayers
%   SSALB_gas      vector, Nlayers
%   NPmom_gas      scalar
%   Pmom_gas       array,  M x Nlayers (M=21, e.g.)
%   DTAUC_wat      vector, Nlayers 
%   SSALB_wat      vector, Nlayers
%   NPmom_wat      scalar
%   Pmom_wat       array,  M2 x Nlayers (M2=30, e.g.)
%
% Outputs:
%   DTAUC_new      vector, Nlayers x 1
%   SSALB_new      vector, Nlayers x 1
%   NPmom_new      scalar
%   Pmom_new       array,  max(M,M2) x Nlayers
%
% - For the vectors, the dimensions may be flipped (doesn't matter)
% - For the matrices (Pmoms) dimensions must be moments x layers
%
% By Steven Neshyba, Aug. 2011
% Comments and formatting added: Penny Rowe, Sept. 8, 2011
% If loops added to check dimensions: PR Sept. 8, 2011
%
% Notes: see sspfuns.py for python version
%

% The variables need to be double precision!
DTAUC_gas = double(DTAUC_gas);  SSALB_gas = double(SSALB_gas); 
NPmom_gas = double(NPmom_gas);  Pmom_gas  = double(Pmom_gas); 
DTAUC_wat = double(DTAUC_wat);  SSALB_wat = double(SSALB_wat); 
NPmom_wat = double(NPmom_wat);  Pmom_wat  = double(Pmom_wat) ;




% Get SSALB as combination of those for gas and water, with layer
[NLYR, Pmom_gas_dim] = size(Pmom_gas) ;
[NLYR, Pmom_wat_dim] = size(Pmom_wat) ;
DTAUC_new = DTAUC_gas + DTAUC_wat ;
SSALB_new = zeros(NLYR,1) ;
NPmom_new = zeros(NLYR,1) ;
for iLYR = 1:NLYR
  % New if/else to catch cases where DTAUC_new (and thus the denom) = 0
  % PR 31 Jan 2012
  if DTAUC_new(iLYR)==0
    SSALB_new(iLYR)=0;
  else
    SSALB_new(iLYR) = ( ...
      SSALB_gas(iLYR)*DTAUC_gas(iLYR) + SSALB_wat(iLYR)*DTAUC_wat(iLYR) )...
      / DTAUC_new(iLYR) ;
  end
end

% Initialize Pmom_new with same size as largest of Pmom_gas and Pmom_wat
if (Pmom_wat_dim > Pmom_gas_dim)
  Pmom_new = zeros(NLYR,Pmom_wat_dim) ;
  Pmom_new_dim = Pmom_wat_dim ;
else
  Pmom_new = zeros(NLYR,Pmom_gas_dim) ;
  Pmom_new_dim = Pmom_gas_dim ;
end

% Get Pmom_new as function of layer, index to moment and NPmom_new(layer)
% Loop over layer
for iLYR = 1:NLYR
  NPmom_new(iLYR) = 0 ;
  denom = DTAUC_new(iLYR)*SSALB_new(iLYR) ;
  %print 'denom: ', denom.shape, denom
  if (denom > 0)
    % Loop over index to Legendre moment
    for iPmom =1:Pmom_new_dim
      %print 'Got a non-zero denom at ', iLYR
      if (iPmom <= NPmom_wat && iPmom <= NPmom_gas)
        % Need to combine Pmoms for both wat and gas
        Pmom_new(iLYR,iPmom) = (...
          SSALB_gas(iLYR)*DTAUC_gas(iLYR)*Pmom_gas(iLYR,iPmom) + ...
          SSALB_wat(iLYR)*DTAUC_wat(iLYR)*Pmom_wat(iLYR,iPmom))/denom ;
        NPmom_new(iLYR) = NPmom_new(iLYR) + 1 ;
        
      elseif (iPmom <= NPmom_wat)
        % Need to combine Pmoms for wat only
        Pmom_new(iLYR,iPmom) = (...
          SSALB_wat(iLYR)*DTAUC_wat(iLYR)*Pmom_wat(iLYR,iPmom))/denom ;
        NPmom_new(iLYR) = NPmom_new(iLYR) + 1 ;
        
      elseif (iPmom <= NPmom_gas)
        % Need to combine Pmoms for gas only
        Pmom_new(iLYR,iPmom) = (...
          SSALB_gas(iLYR)*DTAUC_gas(iLYR)*Pmom_gas(iLYR,iPmom))/denom ;
        NPmom_new(iLYR) = NPmom_new(iLYR) + 1 ;
        
      else
        break
        %print 'ran out of moments for layer ', iLYR
      end
    end
  end
end
NPmom_new = max(NPmom_new) ;




