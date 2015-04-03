function [qext_i,w0_i,NPmom_i,Pmom_i] = sspinterp3( ...
    reff_mesh,wnum_mesh,qext_mesh,w0_mesh,NPmom_mesh,Pmom_mesh,reff_i,wnum_i)
%size(reff_mesh)
%size(w0_mesh)
%size(Pmom_mesh)
%size(reff_i)
%size(wnum_i)
%size(NPmom_mesh)
debug = 0;

if(debug == 1)
    disp('doing nothing')
    qext_i=[];
    w0_i=[];
    NPmom_i=[];
    Pmom_i=[];
else


reff_vec = reff_mesh(:,1)'; %size(reff_vec), reff_vec(1:5)
wnum_vec = wnum_mesh(1,:);  %size(wnum_vec), wnum_vec(1:5)


% First call gets the weights and indices
[qext_i, s, t, onemt, ndx, nrows] = interp2_modified(...
  wnum_vec,reff_vec,qext_mesh,wnum_i, reff_i);


% Subsequent calls reuse the weights and indices
%w0_i = interp2(wnum_vec,reff_vec,w0_mesh,  wnum_i, reff_i); 
% This was a check
w0_i = fagain(w0_mesh, s, t, onemt, ndx, nrows);

% Even get an interpolated number of moments!
NPmom_i = round(fagain(NPmom_mesh, s, t, onemt, ndx, nrows));
NPmom_max = max(NPmom_i);

[M,dum3] = size(reff_i);
Pmom_i = zeros(M,NPmom_max);
temp = [];
for j=1:NPmom_max
    %Pmom_mesh_j = squeeze(Pmom_mesh(:,:,j)); 
    % Not actually necessary to squeeze
    %temp = interp2(wnum_vec,reff_vec,Pmom_mesh_j,  wnum_i, reff_i) 
    % This was a check
    %temp = interp2(wnum_vec,reff_vec,Pmom_mesh(:,:,j),  wnum_i, reff_i) 
    % This was a check
    %temp = fagain(Pmom_mesh_j, s, t, onemt, ndx, nrows)
    temp = fagain(Pmom_mesh(:,:,j), s, t, onemt, ndx, nrows);
    Pmom_i(:,j) = temp;
end
end

function F = fagain(arg3, s, t, onemt, ndx, nrows)
F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
         ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;


     
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       