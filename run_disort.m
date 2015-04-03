
% % % % % % % % % % % % % % % % % % % % % %
function [raddn,radup, rfldn, flup,izm] = run_disort(nus, clrod, ...
  reff_wat, cldODvis_wat, cldlyr_wat,...
  reff_ice, cldODvis_ice, cldlyr_ice,...
  sspWat,   iTemp_wat, wTemp_wat, ...
  sspIce,   iTemp_ice, wTemp_ice, ...
  nlyr, temper, nstr, lamber, phi0, umu0, umu, albedo, ...
  fbeam_nu, delv_nu, procNoStr, debugflag,dirdisort)
% % % % % % % % % % % % % % % % % % % % % %


%
% Written by Steven Neshyba and Penny M. Rowe
%
% Notes: 
%   Inputs temper and clrod must be from the TOA down! 
%


% Lengths of variables & other conveniences
nnus = length(nus);
Nwats = int32(length(sspWat));
Nices = int32(length(sspIce));


% Stuff for explicit variables
NcldLYR_wat = int32(size(iTemp_wat,1));
NcldLYR_ice = int32(size(iTemp_ice,1));


% fixed disort parameters
plank  = true;           % true => thermal emission from layers.
usrtau = true;           % false => return results at all layer boundaries
%                        % true  => return results at specified taus
usrang = true;           % true => specify observation zenith angle
                         % obs zenith angles, neg = dnwllng, in order!
numu   = length(umu);    % number of observation zenith angles
nphi   = 1;              % there will be one observation azimuthal angle
phi    = 0;              % observation azimuthal angle
ibcnd  = 0;              % boundary condition flag, 0=> thermal emis
onlyfl = false;          % 0 => return azimuthally averaged intensities
%                        % (not just fluxes)
fisot  = 0;              % top-bound isotropic illum.: w/m2 if planck=true
btemp  = temper(nlyr+1); % bottom temperature (1 more temp than layers)
ttemp  = temper(0+1);    % top temperature (irrelevant if we set top
%                        % emissivity, temis, equal to zero)
temis  = 0;              % top emissivity
%                        % 0 at TOA for upwelling.  
%                        % sum(dtauc_tot) for downwelling (umu<0), at srfc
%accur	   = 0.0001;     % should be between 0 and 0.01
%header     = 0;




iangdn = find(umu<0) ; Ndn = length(iangdn);
iangup = find(umu>0) ; Nup = length(iangup);



% ... qc: make sure inputs are ok
if umu0>1; error('umu0 out of range (greater than 1'); end
%if find(umu<-1) | find(umu>1); warning('umu out of range?'); end
if umu0>0 && length(fbeam_nu) ~= nnus
  error('fbeam_nu must have same nu vector as ods.');
end
if length(albedo) ~= nnus
  error('albedo must have same nu vector as ods.');
end



% the cloud model (water)
cldytau_vis_wat = zeros(Nwats,nlyr);
for icldLYR_wat = 1:NcldLYR_wat
  iallwats = find(wTemp_wat(icldLYR_wat,:)>0) ;  % w>0
  for i = 1:length(iallwats)
    iwats = iTemp_wat(icldLYR_wat,iallwats(i));
    cldytau_vis_wat(iwats,cldlyr_wat(icldLYR_wat)) = ...
      cldODvis_wat(icldLYR_wat) * wTemp_wat(icldLYR_wat,iallwats(i));
  end
end
if (debugflag)
  disp('Debugging: optical depths for water');
  %[0 ssp_watTemps; temper(1:end-1) cldytau_vis_wat']
  disp([0  temper(1:end-1) cldytau_vis_wat']);
end
cldytau_vis_wat_usingit = max(cldytau_vis_wat,[],2)';
iwat_usingit = find(cldytau_vis_wat_usingit);

        
% the cloud model (ice)
cldytau_vis_ice = zeros(Nices,nlyr);
for icldLYR_ice = 1:NcldLYR_ice
  iallices = find(wTemp_ice(icldLYR_ice,:)>0) ;
  for i = 1:length(iallices) 
    iices = iTemp_ice(icldLYR_ice,iallices(i));
    cldytau_vis_ice(iices,cldlyr_ice(icldLYR_ice)) = ...
      cldODvis_ice(icldLYR_ice) * wTemp_ice(icldLYR_ice,iallices(i));
  end
end
if (debugflag)
  disp('Debugging: optical depths for ice');
  %[0 ssp_iceTemps; temper(1:end-1) cldytau_vis_ice']
  disp([0 temper(1:end-1) cldytau_vis_ice']);
end
cldytau_vis_ice_usingit = max(cldytau_vis_ice,[],2)';
iice_usingit = find(cldytau_vis_ice_usingit);



% % If there is no cloud, add a very thin water cloud so disort won't bomb
% if sum(cldODvis_ice) + sum(cldODvis_wat) == 0
%   warning('No cloud!');
%   %cldODvis_ice(1,1) = 0.00001;
% end


% set up the reff arrays appropriate for this microwindow
reff_microwindow_wat = ones(nnus,1)*reff_wat ;
reff_microwindow_ice = ones(nnus,1)*reff_ice ;


% interpolate water optical properties for this microwindow
% ... output are all vectors of size nnus x 1
% ... except cldypmom_wat, which is nnus x max(cldynpmom_wat)

        
        
% interpolate liq optical properties for microwindow (only if using it)
itot = 0 ;
for iwats = iwat_usingit
  itot = itot+1 ;
    [cldyqext_wat,cldyw0_wat,cldynpmom_wat,cldypmom_wat] = sspinterp3(...
      sspWat(iwats).reffM, ...
      sspWat(iwats).wnumM, ...
      sspWat(iwats).qextM, ...
      sspWat(iwats).w0M,...
      sspWat(iwats).NpmomM, ...
      sspWat(iwats).PmomM, ...
      reff_microwindow_wat,nus);
    total_wat(itot).cldyqext_wat = cldyqext_wat;
    total_wat(itot).cldyw0_wat = cldyw0_wat;
    total_wat(itot).cldynpmom_wat = cldynpmom_wat;
    total_wat(itot).cldypmom_wat = cldypmom_wat;
end

% interpolate ice optical properties for microwindow (only if using it)
itot = 0 ;
for iices = iice_usingit
  itot = itot+1;
  [cldyqext_ice,cldyw0_ice,cldynpmom_ice,cldypmom_ice] = sspinterp3(...
    sspIce(iices).reffM, ...
    sspIce(iices).wnumM, ...
    sspIce(iices).qextM, ...
    sspIce(iices).w0M,...
    sspIce(iices).NpmomM, ...
    sspIce(iices).PmomM, ...
    reff_microwindow_ice,nus);
  total_ice(itot).cldyqext_ice = cldyqext_ice;
  total_ice(itot).cldyw0_ice = cldyw0_ice;
  total_ice(itot).cldynpmom_ice = cldynpmom_ice;
  total_ice(itot).cldypmom_ice = cldypmom_ice;
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% ... initialize the radiance
if Ndn>0; raddn = zeros(nnus,Ndn) ; else; raddn=[];
if Nup>0; radup = zeros(nnus,Nup) ; else; radup=[];
rfldn = zeros(nnus,1) ;
flup  = zeros(nnus,1) ;


% solar contribution only if umu0>0
if umu0<=0
  % sun is below the horizon, so turn it off
  fbeam=0;
  umu0=0;
end


% loop over all the wavenumbers
izm = zeros(length(nus),1);   % added Sept. 9, 2012, PMR
for inu = 1:length(nus)
  
  pause(0)
  
  % width of microwindow for disort input (fluxes) pmr 2 feb. 2012
  delv = delv_nu(inu) ;
  
  % wavenumber for this iteration
  wvnumlo = nus(inu)-delv/2 ;  % uses delv pmr 2 feb. 2012
  wvnumhi = nus(inu)+delv/2 ;
  
  % solar contribution: only if umu0>0
  if umu0>0
    fbeam = fbeam_nu(inu) ;
  end
  
  % clear sky optical depth contribution for this frequency
  dtauc_tot = clrod(:,inu)' ;          % clear sky optical depth (nnus x 1)
  ssalb_tot = zeros(nlyr,1) ;          % 0s (1 x nlyr)
  npmom_tot = 0 ;
  pmom_tot  = zeros(nlyr,1) ;          % 0s (nmom_tot+1,nlyr)
  
  % cloud water contribution for this frequency (only if we're using it)
  itot = 0 ;
  for iwats = iwat_usingit
    itot = itot+1;
    dtauc_wat_inu = total_wat(itot).cldyqext_wat(inu)*...
      cldytau_vis_wat(iwats,:)/2;                                % 1 x nlyr
    ssalb_wat_inu = total_wat(itot).cldyw0_wat(inu)*ones(nlyr,1);% nlyr x 1
    npmom_wat_inu = total_wat(itot).cldynpmom_wat(inu);          % 1 x 1
    pmom_wat_inu  = ones(nlyr,1)*...
      total_wat(itot).cldypmom_wat(inu,1:npmom_wat_inu);         % ?
    [dtauc_tot,ssalb_tot,npmom_tot,pmom_tot] = combine2(...
      dtauc_tot,ssalb_tot,npmom_tot,pmom_tot,...
      dtauc_wat_inu,ssalb_wat_inu,npmom_wat_inu,pmom_wat_inu);
  end
  
                                                 
  % cloud ice contribution for this frequency (only if we're using it)
  itot = 0 ;
  for iices = iice_usingit
    itot = itot+1 ;
    dtauc_ice_inu = total_ice(itot).cldyqext_ice(inu)*cldytau_vis_ice(iices,:)/2; % 1 x nlyr
    ssalb_ice_inu = total_ice(itot).cldyw0_ice(inu)*ones(nlyr,1);  % nlyr x 1
    npmom_ice_inu = total_ice(itot).cldynpmom_ice(inu); % 1 x 1
    pmom_ice_inu  = ones(nlyr,1)*total_ice(itot).cldypmom_ice(inu,1:npmom_ice_inu);
    [dtauc_tot,ssalb_tot,npmom_tot,pmom_tot] = combine2(...
      dtauc_tot,ssalb_tot,npmom_tot,pmom_tot,...
      dtauc_ice_inu,ssalb_ice_inu,npmom_ice_inu,pmom_ice_inu);
  end
  
  % ... Cut out layers at top where od<5e-7 or 1e-5, otherwise round-off error
  %     can become significant
  i2 = length(dtauc_tot);
  izeds = find(dtauc_tot<1e-5) ;
  if isempty(izeds)
    i1 = 1 ;
  elseif izeds(1)~=1 || any(diff(izeds)~=1)
    %warning('Redesign your layers so od gets smaller going up.');
    %izeds
    i1 = max(izeds)+1 ;
  else
    i1 = max(izeds)+1 ;
  end
  izm(inu)=i1;
  
  % And get ssalb_totNew for the new atmospheric profile
  ssalb_totNew = ssalb_tot(i1:i2)  ;
  
  % % % % % % % % % % % % % % % % % % % % %
  % ... extra code needed to run from matlab
  
  % Locate any layers with nonzero ssalb
  cldlyr     = int32(find(ssalb_tot));
  cldlyrNew  = int32(find(ssalb_totNew));
  Ncldlyr    = int32(length(cldlyr));
  
  
  % cumulative optical depth of atmosphere
  % for downwelling, utau is the cumulative optical depth of atmosphere
  % while for upwelling leave utau set to 0
  
  % Added "single"!  PMR Sept. 6, 2012
  utau(1) = sum(single(dtauc_tot(i1:i2))) ;
  utau(2) = 0 ;
  % if isempty(cldlyrNew)
  %   utau(1) = sum(dtauc_tot(i1:cldlyr_wat)) ;
  %   utau(2) = sum(dtauc_tot(i1:cldlyr_wat-1)) ;
  % else 
  %   utau(1) = sum(dtauc_tot(i1:cldlyrNew)) ;
  %   utau(2) = sum(dtauc_tot(i1:cldlyrNew-1)) ;
  % end
  
  %utau = sum(dtauc_tot) ;
  %utau = 0 ;
  ntau   = length(utau);    % number of specified taus (output heights)

  
  % Open the namelist file
  fid=fopen([dirdisort 'disortinput.nml'],'w');
  
  % write disort namelist
  fprintf(fid,'%s\n','&disortinput');
  write_namelist_int(fid,'nstr',  nstr);
  write_namelist_int(fid,'nlyr',  i2-i1+1);
  write_namelist_int(fid,'ntau',  ntau);
  write_namelist_vec(fid,'utau',  utau);
  write_namelist_int(fid,'nphi',  nphi);
  write_namelist_vec(fid,'phi',   phi);
  write_namelist_int(fid,'ibcnd', ibcnd);
  write_namelist_flt(fid,'umu0',  umu0);
  write_namelist_flt(fid,'phi0',  phi0);
  write_namelist_flt(fid,'albedo',albedo(inu));
  write_namelist_vec(fid,'ssalb', ssalb_totNew);
  write_namelist_flt(fid,'fbeam', fbeam);
  write_namelist_flt(fid,'fisot', fisot);
  write_namelist_vec(fid,'dtauc', dtauc_tot(i1:i2));
  write_namelist_flt(fid,'wvnmlo',wvnumlo);
  write_namelist_flt(fid,'wvnmhi',wvnumhi);
  write_namelist_boo(fid,'usrtau',usrtau);
  write_namelist_int(fid,'numu',  numu);
  write_namelist_vec(fid,'umu',   umu);
  write_namelist_boo(fid,'usrang',usrang);
  write_namelist_boo(fid,'lamber',lamber);
  write_namelist_flt(fid,'temis', temis);
  write_namelist_boo(fid,'plank', plank);
  write_namelist_boo(fid,'onlyfl',onlyfl);
  write_namelist_boo(fid,'debug', 0);
  write_namelist_vec(fid,'temper',temper(i1:i2+1));
  write_namelist_vec(fid,'ttemp', temper(i1));
  write_namelist_vec(fid,'btemp', btemp);
  write_namelist_int(fid,'nmom',  npmom_tot-1);
  write_namelist_int(fid,'Ncldlyr',Ncldlyr);
  outfile = 'disort_out.txt';
  write_namelist_str(fid,'outfile',outfile);
  fprintf(fid,'%s\n\n','/');
  
  % Tack on moments afterward
  for icldlyr = 1:Ncldlyr
    
    % Failsafe check: If pmom_cld(1) >= 1.00000006, disort will complain
    if pmom_tot(icldlyr,1)>=1.00000006; error('pmom(1) is too big!'); end
    
    
    % Write out the moments for this layer
    fprintf(fid,'%s\n','&PMOMINPUT');
    try
      write_namelist_int(fid,'CLDLYR',cldlyrNew(icldlyr));
    catch
      dum=temper(i1:i2+1); disp([temper(cldlyr) dum(cldlyrNew)]);
      error('Could not write new cloud layer to namelist.');
    end
    write_namelist_vec(fid,'PMOM_CLD',pmom_tot(cldlyr(icldlyr),:));
    fprintf(fid,'%s\n\n','/');
    
  end
  
  % Close the namelist file
  fclose(fid);
  %type ([dirdisort 'disortinput.nml'])
  
  % % % % % % % % % % % % % % % % % % % % %
  
  trfldn = nan ; tflup=nan; % this logic b/c disort doesn't always converge
  while isnan(trfldn) || isnan(tflup)
    pause(0)
    % run disort
    % note: if you get a warning that the temperature change may be too
    % large for good accuracy, fix the problem, because it is true!)
    cd(dirdisort);         % running disort takes about 0.08 s
                                                 
    % .. Run DISORT using the fortran code for Matlab: disort_driver_mat
    %    the call to system is different for Octave and Matlab
    if (is_octave)
      % Octave
      system('./disort_driver_mat');
    else
      % matlab
      eval('!./disort_driver_mat');
    end
                                                 
    
    % ... Get results
    result = load(outfile) ;
    trfldn = result(1);
    tflup = result(2);
    %if find(uudn<0)
    %error('uu less than zero');
    %end
  end
  uudn = result(3:3+numu-1);
  uuup = result(3+numu:3+numu+numu-1);
  for k=1:Ndn
    raddn(inu,k) = uudn(iangdn(k)) / (wvnumhi - wvnumlo) ;
  end
  for k=1:Nup
    radup(inu,k) = uuup(iangup(k)) / (wvnumhi - wvnumlo) ;
  end
  rfldn(inu) = result(1) / (wvnumhi - wvnumlo) ;
  flup(inu)  = result(2) / (wvnumhi - wvnumlo) ;
  %disp (['done w/ ', num2str(inu) ' of ' num2str(length(nus))])
end
%unix(['rm -f ',outfile,' disortinput.nml']);
end

