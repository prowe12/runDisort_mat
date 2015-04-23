function [msg,x,y,z,out5,out6] = xyzchk(arg1,arg2,arg3,arg4,arg5,arg6)
%XYZCHK Check arguments to 3-D data routines.
%	[MSG,X,Y,Z,C] = XYZCHK(Z), or
%	[MSG,X,Y,Z,C] = XYZCHK(Z,C), or
%	[MSG,X,Y,Z,C] = XYZCHK(X,Y,Z), or
%	[MSG,X,Y,Z,C] = XYZCHK(X,Y,Z,C), or
%	[MSG,X,Y,Z,XI,YI] = XYZCHK(X,Y,Z,XI,YI) checks the input aguments
%	and returns either an error message in MSG or valid X,Y,Z (and
%	XI,YI) data.
%
%	[MSG,X,Y,Z,XI,YI] = XYZCHK(X,Y,Z,XI,YI,'map') maps data 
%	coordinates (X,Y) and (XI,YI) to an equally spaced domain
%	when data domain is not equally spaced.  Note that (X,Y) 
%	are assumed to be plaid in this case.

%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.10 $
%
%   Modified for Octave!


warning('off', 'Octave:possible-matlab-short-circuit-operator');

error(nargchk(1,6,nargin));

msg = [];

if nargin==1, % xyzchk(z)
  z = arg1;
  if isstr(z),
    msg = 'Input arguments must be numeric.';
    return
  end
  [m,n] = size(z);
  [x,y] = meshgrid(1:n,1:m);
  out5 = z; % Default color matrix
  return

elseif nargin==2, % xyzchk(z,c)
  z = arg1; c = arg2;
  [m,n] = size(z);
  [x,y] = meshgrid(1:n,1:m);
  if any(size(z)~=size(c)),
    msg = 'Matrix C must be the same size as Z.';
    return
  end
  out5 = c;
  return

elseif nargin>=3, % xyzchk(x,y,z,...)
  x = arg1; y = arg2; z = arg3;
  [m,n] = size(z);
  if min(size(z))>1, % z is a matrix
    % Convert x,y to row and column matrices if necessary.
    % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
    if min(size(x))==1 && min(size(y))==1,
      [x,y] = meshgrid(x,y);
      % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
      %if size(x,2)~=n | size(y,1)~=m,
      if size(x,2)~=n || size(y,1)~=m,
        msg = 'The lengths of X and Y must match the size of Z.';
        return
      end
    % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
    %elseif min(size(x))==1 | min(size(y))==1,
    elseif min(size(x))==1 || min(size(y))==1,
      msg = 'X and Y must both be vectors or both be matrices.';
      return
    else
      % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
      %if any(size(x)~=size(z)) | any(size(y)~=size(z)),
      if any(size(x)~=size(z)) || any(size(y)~=size(z)),
        msg = 'Matrices X and Y must be the same size as Z.';
        return
      end
    end
  else % z is a vector
    % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
    %if min(size(x))~=1 | min(size(y))~=1,
    if min(size(x))~=1 || min(size(y))~=1,
      msg = 'X and Y must be vectors when Z is.';
      return
    % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
    %elseif length(x)~=length(z) | length(y)~=length(z),
    elseif length(x)~=length(z) || length(y)~=length(z),
         msg = 'X and Y must be same length as Z.';
      return
    end
  end
end

if nargin==4, % xyzchk(x,y,z,c)
  c = arg4;
  if any(size(z)~=size(c)),
    msg = 'Matrix C must be the same size as Z.';
    return
  end
  out5 = c;
  return
end

if nargin>4, % xyzchk(x,y,z,xi,yi,...)
  xi = arg4; yi = arg5;
end

% Check for non-equally spaced data.  If so, map (x,y) and
% (xi,yi) to matrix (row,col) coordinate system.
if (nargin==6)
  if isstr(arg6), % xyzchk(x,y,z,xi,yi,'map')
    xx = x(1,:); yy = y(:,1);
    dx = diff(xx); dy = diff(yy);
    % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
    if (max(abs(diff(dx))) > eps*max(xx)) || ...
       (max(abs(diff(dy))) > eps*max(yy)),
      if any(dx < 0), % Flip orientation of data so x is increasing.
        x(:) = fliplr(x); y(:) = fliplr(y); z(:) = fliplr(z);
        xx(:) = fliplr(xx); dx(:) = -fliplr(dx);
      end
      if any(dy < 0), % Flip orientation of data so y is increasing.
        x(:) = flipud(x); y(:) = flipud(y); z(:) = flipud(z);
        yy(:) = flipud(yy); dy(:) = -flipud(dy);
      end

      % .. converted to short-circuit operator for Octave, PMR, 2015/23/01
      if any(dx<=0) || any(dy<=0),
        msg = 'X and Y must be monotonic.';
        return
      end

      % Map values in xi to values in ui via linear interpolation.
      ui = xx(1)*ones(size(xi)); % Initialize ui with out of range value
      d = (xi < xx(1));
      for i=1:length(xx)-1,
        e = (xi <= xx(i+1)) & ~d;
        f = find(e);
        if ~isempty(f), ui(f) = i + (xi(f)-xx(i))/dx(i); end
        d = d | e;
      end

      % Map values in yi to values in vi via linear interpolation.
      vi = yy(1)*ones(size(yi)); % Initialize vi with out of range value
      d = (yi < yy(1));
      for i=1:length(yy)-1,
        e = (yi <= yy(i+1)) & ~d;
        f = find(e);
        if ~isempty(f), vi(f) = i + (yi(f)-yy(i))/dy(i); end	
        d = d | e;
      end

      [x,y] = meshgrid(1:size(x,2),1:size(y,1));
      xi = ui; yi = vi;
    end
  end
end


if nargin>4, % Check interpolation arguments
  % If x is a row and y is a column (or vice versa) then build xi,yi matrices.
  % .. converted to short-circuit operator for Octave, PMR, 2015/03/01
  if min(size(xi))==1 && min(size(yi))==1 && any(size(xi)~=size(yi))
    [xi,yi] = meshgrid(xi,yi);
  elseif any(size(xi)~= size(yi)), % Also create matrix if sizes differ
    [xi,yi] = meshgrid(xi,yi);
  end
  out5 = xi; out6 = yi;
end