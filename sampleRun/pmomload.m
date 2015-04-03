function [Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,...
  qext_mesh,wnum_list,reff_list] = pmomload(pmomfile)
%
%
% function [Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,...
%          qext_mesh,wnum_list,reff_list] = pmomload(pmomfile)
%
% Purpose: Load pmoment data from file
%
% by Penny M. Rowe and Steven Neshyba
%
% Last update: 2015/03/01
%   Added optional line to enable netcdf in octave
%

% ... If working in Octave, need to import netcdf package
if (is_octave)
  import_netcdf
end

nc = netcdf.open(pmomfile,'NC_NOWRITE') ;
    
    
% Flip all of these to be comparable to python code
% (matlab get the transpose compared to python)
Npmomarray = netcdf.getVar(nc,netcdf.inqVarID(nc,'Npmomarray'))';
pmomarray0 = netcdf.getVar(nc,netcdf.inqVarID(nc,'pmomarray'));
wnum_mesh  = netcdf.getVar(nc,netcdf.inqVarID(nc,'wnum_mesh'))';
reff_mesh  = netcdf.getVar(nc,netcdf.inqVarID(nc,'reff_mesh'))';
w0_mesh    = netcdf.getVar(nc,netcdf.inqVarID(nc,'w0_mesh'))';
qext_mesh  = netcdf.getVar(nc,netcdf.inqVarID(nc,'qext_mesh'))';
reff_list  = netcdf.getVar(nc,netcdf.inqVarID(nc,'reff_list'))';
wnum_list  = netcdf.getVar(nc,netcdf.inqVarID(nc,'wnum_list'))';

% Fix dimensions of pmomarray, because matlab get different dimensions
% than python.  I found the correct dims by trial and error.
[d1,d2,d3] = size(pmomarray0);
pmomarray = zeros(d3,d2,d1);
for i=1:d2
  pmomarray(:,i,:) = squeeze(pmomarray0(:,i,:))';
end

% ... Make double, not int32, for matlab
Npmomarray = double(Npmomarray) ;

    
% Npmomarray = nc.variables['Npmomarray'][:]; Npmomarray = Npmomarray.astype('int32')
% pmomarray  = nc.variables['pmomarray'][:]; pmomarray = pmomarray.astype('float64')
% wnum_mesh  = nc.variables['wnum_mesh'][:]; wnum_mesh = wnum_mesh.astype('float64')
% reff_mesh  = nc.variables['reff_mesh'][:]; reff_mesh = reff_mesh.astype('float64')
% w0_mesh    = nc.variables['w0_mesh'][:]; w0_mesh = w0_mesh.astype('float64')
% qext_mesh  = nc.variables['qext_mesh'][:]; qext_mesh = qext_mesh.astype('float64')
% reff_list  = nc.variables['reff_list'][:]; reff_list = reff_list.astype('float64')
% wnum_list  = nc.variables['wnum_list'][:]; wnum_list = wnum_list.astype('float64')

netcdf.close(nc)

    
    
