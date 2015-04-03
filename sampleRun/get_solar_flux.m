
function [ssf_flux] = get_solar_flux(ssf_wnum,ssf_spectra,distance,...
  wnum1,wnum2)

ssf_spectra1 = interp1(ssf_wnum,ssf_spectra,wnum1);
ssf_spectra2 = interp1(ssf_wnum,ssf_spectra,wnum2);

% Find the two indices that bound the passband specified by wnum1 and wnum2
is = find(ssf_wnum >= wnum1, 1);        % index to starting nu (before)
ie = find(ssf_wnum <= wnum2, 1,'last'); % index to ending nu (after)

ssf_flux = trapz([wnum1 ssf_wnum(is:ie) wnum2],...
  [ssf_spectra1 ssf_spectra(is:ie) ssf_spectra2]) ;

ssf_flux = ssf_flux/ distance^2 ;

