function [Iout,wout] = interp1_weights2(y,yi_list,method)

% Purpose: Returns indices and weights
%
% By S. Neshyba and P. Rowe
%

% Sort the y's in increasing order
[ysorted, iysorted] = sort(y);
Ny = length(y);                        % Number of ys
I = int32([]);
w = [];
N_yi_list = length(yi_list);           % Number of interpolated ys (yis)

for iyi = 1:N_yi_list
  
  yi = yi_list(iyi);                   % Current interpolated y (yi)
  
  if (yi >= ysorted(Ny))               % yi is bigger than biggest y
    I1 = Ny; w1 = 1;
    I2 = Ny; w2 = 0;   % PMR 7 Mar 2012
  elseif (yi <= ysorted(1))            % yi is smaller than smallest y
    I1 = 1; w1 = 1;
    I2 = 1; w2 = 0;    % PMR 7 Mar 2012
  else
    Ifloat = interp1(ysorted, [1:Ny], yi);
    I1 = floor(Ifloat);
    I2 = ceil(Ifloat);
    if(I2==I1)
      I2=I1+1;
    end
    w1 = (I2-Ifloat)/(I2-I1);
    w2 = 1-w1;
  end
  I = [I;[iysorted(I1) iysorted(I2)]];
  w = [w;[w1 w2]];
  
end

if (nargin > 2)
  if (strcmp(method,'nearest'))
    [wnearest, iwnearest] = max(w'); wnearest = wnearest';
    Inearest = int32(zeros(N_yi_list,1));
    for iyi = 1:N_yi_list
      Inearest(iyi,1) = I(iyi,iwnearest(iyi));
    end
    Iout = Inearest;
    wout = ones(size(Inearest));
    if length(Iout)>1; error('Check code: length(Iout)>1.'); end
    elseif (strcmp(method,'interpolate'))
      Iout = I;
      wout = w;
    else
      error(['Method should be ''nearest'' or ''interpolate,'' ' ...
        'but is ' method '.']);
      Iout = [];
      wout = [];
    end
  end

end

