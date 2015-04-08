
% sample_run.m
%
% Purpose: A sample code for running DISORT, using dummy wavenumbers
%   and optical depths given here.
%
% Inputs (see below)
%
% Outputs (will change to allow user to specify viewing angle)
%   raddn: Downwelling radiances for 180 degrees
%          (180 => looking straight up)
%   radup: Upwelling radiance (not included here)
%          (0 => looking straight down)
%   rfldn: Downward flux
%   flup:  Upward flux
%   izm:   index to top layer used (cuts off where od<1e-5), where 1=>top
%
% Penny M. Rowe (prowe@harbornet.com)
% July 18, 2004
%
% Note: can be run on either Matlab or Octave
%



% % % % % % % %   INPUTS (replace as desired)   % % % % % % % % % % %

% ... Directories
dirdisort = '../runDisort_mat/';       % directory where output will go (for now, here)
sspdir    = '../runDisort_mat/inputs/';

% ... Path needed for Octave
if (is_octave)
  addpath '../runDisort_mat/octCode/'
end




% ... Legendre moment files for water and ice, based on single scattering
%     parameters.  
%
%     For water, these are based on temperature-dependent measurements 
%     of the refractive index:
%     Downing and Williams, 1975 for 300 K
%     Zasetsky et al, 2005 and Wagner et al., 2005 for 240 to 263 K
%
%     For ice, these are based on Warren et al., 1984 and assume spheres
%
%     We hope to add more options for both water and ice soon. 
%
%     For water, we interpolate to the cloud temperature, as indicated
%     by variables below, while for ice we assume temperature-independent
%     (we hope to add habits and temperature-dependent ice eventually as well).
%
pmomfiles_wat = {...
  [sspdir '/liqWater_T240_S331_pmom.nc'],...
  [sspdir '/liqWater_T253_S331_pmom.nc'],...
  [sspdir '/liqWater_T263_S331_pmom.nc'],...
  [sspdir '/liqWater_T300_S100_pmom.nc']};
pmomfiles_ice = {[sspdir '/iceWater_spheres_S100_pmom.nc']};
watTdependence = 'interpolate';  % char string
iceTdependence = 'none';  % char string
ssp_watTemps = [240 253 263 273]; % temp of ssp file, in K
ssp_iceTemps = 273; % temp of spp file, in K


% ... Cloud, atmospheric, and surface parameters for DISORT
dateno = datenum(2014,7,18,12,0,0);
location.latitude = 47.2;     % North of equator is positive
location.longitude = -122.5;  % East of Greenwhich is positive
location.altitude = 243 ;     % meters

% Scene viewing angle
sceneAngle = 180;             % degrees, 180=>looking up (dwnwllng), 0=>down

% Functions of Wavenumber (vectors)
nus    = [850; 950];
albedo = [0.02; 0.02];

% ... Functions of Wavenumber and Layers (matrix).
%     The clear-sky optical depths (gases only) are a function of both
%     frequency and layer (TOA down to surface), so this should be
%     a matrix.  It must have the correct dimensions
clrod   = [...
  .2 .4
  .3 .3
  .2 .5];           % layers x wavenumbers

% ... Function of layer boundaries, so length is # layers + 1
%     Layer temperatures, TOA down to surface
temper  = [273 263 253 243];      % Kelvin, vector


% Integers
nlyr    = size(clrod,1);        % Number of layers (in clrod)
nstr    = 16 ;                  % Number of streams
lamber  = 1 ;                   % Lambertian surface flag
solarDistance = 1 ;

% Cloud properties, one per cloud layer
reff_wat = 10;                  % effective radius (liquid) water
reff_ice = 50;                  % effective radius ice (water)
cldODvis_wat = 1 ;              % optical depth in visible of water
cldODvis_ice = .5;              % optical depth in visible of ice
cldlyr_wat = 1;                 % layer #s containing water cloud, from top
cldlyr_ice = 1;                 % layer #s containing ice cloud, from top
cldTemp    = mean(temper(3:4)); % Temperature of cloud


% Other flags
procNoStr = 1;        % For multiple processors 
debugflag = 0;        % If 1, returns debugging info




% % % % % %    Set up for run    % % % % % % % % % % % %

% spectral spacing (best to leave as 1)
delv   = ones(size(nus));       


% ... Get the Solar zenith angle (SZA)
% Get the location, for getting SZA
sun  = sun_position(dateno, location);
sza  = sun.zenith ;
saa  = sun.azimuth ;
umu0 = cos(sza*pi/180);

% ... Cosine of scene viewing angle and SZA
%     Note: you could have more than one, and/or both
%     up and downwelling, e.g.
%     umu    = [-1 -.5 -0.01 0.01 0.5 1];  
umu = cos(sceneAngle*pi/180) ;  % -1=>looking up (dwnwllng), 1=>down



% ... Solar spectrum
nuSun = floor(min(nus))-10:1:ceil(max(nus))+10 ;
radSun= getsolarbeam_IR(nuSun); % Incident beam intensity, W m^-2 cm

% ... Get solar input, fbeam:
fbeam = 0*nus ;
for iinu = 1:length(nus)
  fbeam(iinu) = get_solar_flux(nuSun,radSun,solarDistance,...
    nus(iinu)-delv(iinu)/2,nus(iinu)+delv(iinu)/2);
end



% ... Get Legendre moments for water and ice
%     Water
for i=1:length(pmomfiles_wat)
  [sspWat(i).NpmomM,sspWat(i).PmomM,sspWat(i).wnumM,sspWat(i).reffM,...
    sspWat(i).w0M,sspWat(i).qextM,sspWat(i).wnum_list,...
    sspWat(i).reff_list] = pmomload(char(pmomfiles_wat(i)));
end

%     Ice
for i=1:length(pmomfiles_ice)
  [sspIce(i).NpmomM,sspIce(i).PmomM,sspIce(i).wnumM,sspIce(i).reffM,...
    sspIce(i).w0M, sspIce(i).qextM,sspIce(i).wnum_list,...
    sspIce(i).reff_list] = pmomload(char(pmomfiles_ice(i)));
end

% Temperature dependence of ssp data: water and ice
[iTempWat,wTempWat] = ...
  interp1_weights2(ssp_watTemps,cldTemp,'interpolate');
[iTempIce,wTempIce] = ...
  interp1_weights2(ssp_iceTemps,cldTemp,'nearest');


% ... Run DISORT
[raddn,radup, rfldn, flup,izm] = run_disort(nus, clrod, ...
  reff_wat, cldODvis_wat, cldlyr_wat,...
  reff_ice, cldODvis_ice, cldlyr_ice,...
  sspWat,   iTempWat, wTempWat, ...
  sspIce,   iTempIce, wTempIce, ...
  nlyr, temper, nstr, lamber, saa, umu0, umu, albedo, ...
  fbeam, delv, procNoStr, debugflag,dirdisort) ; 

plot(nus,raddn,'.-')
xlabel('wavenumber (cm^-^1)')
ylabel('Radiance [(mW/(m^2 sr cm^-^1)]')
legend([num2str(sceneAngle) '^o'])

fprintf('\n Success! \n')


