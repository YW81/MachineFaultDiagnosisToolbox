function varargout = mySVD(mat, varargin)
%--------------------------------------------------------------------------
% Syntax:       [U S V] = mySVD(mat);
%               [U S V] = mySVD(mat,'full');
%               [U S V] = mySVD(mat,'compact');
%               s = mySVD(mat);
%               s = mySVD(mat,'full');
%               s = mySVD(mat,'compact');
%
% Inputs:       mat is an arbitrary M x N matrix
%               mode can be {'full','compact'}
%
% Outputs:      When mode == 'full':
%               U is an M x M orthogonal matrix, S is an M x N diagonal
%               matrix of singular values, and V is an N x N orthoganl
%               matrix. Alternatively, s is a vector of length min(M,N) of
%               singular values.
%
%               When mode == 'compact':
%               U is an M x R matrix, S is an R x R diagonal matrix of
%               singular values, V is an N x R matrix. Alternatively, s is
%               a vector of length R of singular values. Here R = rank(mat)
%
%               In either case, U, S, and V satisfy the following relation:
%               mat = U * S * V';
%              
% Description:  This function computes the singular value decomposition
%               (SVD) of input mat. When mode == 'compact', this function
%               returns the compact SVD containing only the R pairs of 
%               vectors corresponding to the R nonzero singular values of
%               input mat. Here R = rank(mat).
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         June 28, 2012
%--------------------------------------------------------------------------
  
  % Maximum number of allowed iterations
  MAX_SVD_ITER = 30;
  
  % Input check
  if nargin > 1
      mode = varargin{1};
  else
      mode = 'full';
  end
  [m n] = size(mat);
  if (n > m)
      [Ut St Vt] = mySVD(mat',mode);
      if (nargout == 3)
          varargout{1} = Vt;
          varargout{2} = St';
          varargout{3} = Ut;
      else
          if (min(size(St)) == 1)
            varargout{1} = St(1,1);
          else
            varargout{1} = diag(St);
          end
      end
      return;
  end

  % Initialize variables
  Ufull = zeros(m);
  Ufull(1:m,1:n) = mat;
  svals = zeros(n,1);
  Vfull = zeros(n);
  vect = zeros(n,1);
  l = 0;
  mn = 0;
  g = 0;
  scale = 0.0;
  norm = 0.0;

  % Householder reduction to bidiagonal form
  for i = 1:n
    l = i + 1;
    vect(i) = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for k = i:m
        scale = scale + abs(Ufull(k, i));
    end
    if (scale ~= 0)
      for k = i:m
        Ufull(k, i) = Ufull(k, i) / scale;
        s = s + Ufull(k, i) * Ufull(k, i);
      end
      f = Ufull(i, i);
      g = -sqrt(s) * sign(f);
      h = f * g - s;
      Ufull(i, i) = f - g;
      for j = l:n
        s = 0.0;
        for k = i:m
            s = s + Ufull(k, i) * Ufull(k, j);
        end
        f = s / h;
        for k = i:m
            Ufull(k, j) = Ufull(k, j) + f * Ufull(k, i);
        end
      end
      for k = i:m
          Ufull(k, i) = Ufull(k, i) * scale;
      end
    end
    svals(i) = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for k = l:n
        scale = scale + abs(Ufull(i, k));
    end
    if (scale ~= 0)
      for k = l:n
        Ufull(i, k) = Ufull(i, k) / scale;
        s = s + Ufull(i, k) * Ufull(i, k);
      end
      f = Ufull(i, l);
      g = -sqrt(s) * sign(f);
      h = f * g - s;
      Ufull(i, l) = f - g;
      for k = l:n
          vect(k) = Ufull(i, k) / h;
      end
      for j = l:m
        s = 0.0;
        for k = l:n
            s = s + Ufull(j, k) * Ufull(i, k);
        end
        for k = l:n
            Ufull(j, k) = Ufull(j, k) + s * vect(k);
        end
      end
      for k = l:n
          Ufull(i, k) = Ufull(i, k) * scale;
      end
    end
    norm = max(norm, (abs(svals(i)) + abs(vect(i))));
  end
  
  % Accumulate right-hand transformations
  for i = n:-1:1
    if (g ~= 0)
      % Double division to avoid possible underflow
      for j = l:n
          Vfull(j, i) = (Ufull(i, j) / Ufull(i, l)) / g;
      end
      for j = l:n
        s = 0.0;
        for k = 1:n
            s = s + Ufull(i, k) * Vfull(k, j);
        end
        for k = l:n
            Vfull(k, j) = Vfull(k, j) + s * Vfull(k, i);
        end
      end
    end
    for j = l:n
        Vfull(i, j) = 0.0;
        Vfull(j, i) = 0.0;
    end
    Vfull(i, i) = 1.0;
    g = vect(i);
    l = i;
  end

  % Accumulate left-hand transformations
  for i = n:-1:1
    l = i + 1;
    g = svals(i);
    for j = l:n
        Ufull(i, j) = 0.0;
    end
    if (g ~= 0)
      g = 1.0 / g;
      for j = l:n
        s = 0.0;
        for k = l:m
            s = s + Ufull(k, i) * Ufull(k, j);
        end
        f = (s / Ufull(i, i)) * g;
        for k = i:m
            Ufull(k, j) = Ufull(k, j) + f * Ufull(k, i);
        end
      end
      for j = i:m
          Ufull(j, i) = Ufull(j, i) * g;
      end
    else
      for j = i:m
          Ufull(j, i) = 0.0;
      end
    end
    Ufull(i, i) = Ufull(i, i) + 1;
  end

  % Diagonalize the bidiagonal form
  for k = n:-1:1 % loop over all singular values
    for iters = 1:MAX_SVD_ITER % loop over allowed iterations
      flag = 1;
      % Test for splitting
      for l = k:-1:1
        mn = l - 1;
        % Note: vect(1) = 0, always
        if ((abs(vect(l)) + norm) == norm)
          flag = 0;
          break;
        end
        if ((abs(svals(mn)) + norm) == norm)
            break;
        end
      end
      if (flag ~= 0)
        % Cancel vect(l), l > 1
        c = 0.0;
        s = 1.0;
        for i = l:k
          f = s * vect(i);
          vect(i) = c * vect(i);
          if ((abs(f) + norm) == norm)
              break;
          end
          g = svals(i);
          h = SafeDistance(f, g);
          svals(i) = h;
          h = 1.0 / h;
          c = g * h;
          s = -f * h;
          for j = 1:m
            y = Ufull(j, mn);
            z = Ufull(j, i);
            Ufull(j, mn) = y * c + z * s;
            Ufull(j, i) = z * c - y * s;
          end
        end
      end
      z = svals(k);
      if (l == k) % We converged!
        % Make singular value nonnegative
        if (z < 0.0)
          svals(k) = -z;
          for j = 1:n
              Vfull(j, k) = -Vfull(j, k);
          end
        end
        break;
      end
      if (iters == MAX_SVD_ITER)
        disp(['mySVD() reached maximum number of iterations: ' num2str(MAX_SVD_ITER)]);
      end

      % Shift from bottom 2 x 2 minor
      x = svals(l);
      mn = k - 1;
      y = svals(mn);
      g = vect(mn);
      h = vect(k);
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = SafeDistance(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + abs(g) * sign(f))) - h)) / x;

      % Perform next QR decomposition
      c = 1.0;
      s = 1.0;
      for j = l:mn
        i = j + 1;
        g = vect(i);
        y = svals(i);
        h = s * g;
        g = c * g;
        z = SafeDistance(f, h);
        vect(j) = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for jj = 1:n
          x = Vfull(jj, j);
          z = Vfull(jj, i);
          Vfull(jj, j) = x * c + z * s;
          Vfull(jj, i) = z * c - x * s;
        end
        z = SafeDistance(f, h);
        svals(j) = z;
        if (z ~= 0) % Note: Rotation can be arbitrary if z = 0
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        end
        f = c * g + s * y;
        x = c * y - s * g;
        for jj = 1:m
          y = Ufull(jj, j);
          z = Ufull(jj, i);
          Ufull(jj, j) = y * c + z * s;
          Ufull(jj, i) = z * c - y * s;
        end
      end
      vect(l) = 0.0;
      vect(k) = f;
      svals(k) = x;
    end
  end

  if strcmpi(mode,'compact')
    % Compute singular value tolerance
    SVD_TOL = max(m,n) * max(svals) * eps(class(mat));

   % Determine rank to working tolerance
    r = 0;
    for i = 1:n
      if (svals(i) > SVD_TOL)
        r = r + 1;
      end
    end
    
    % Sort ascending
    [svals inds] = sort(svals);
    
    % Return compact SVD
    if (nargout ~= 3)
        varargout{1} = svals(r:-1:1);
    else
        % Populate compact r-dimensional entries
        U = zeros(m, r);
        S = zeros(r);
        V = zeros(n, r);
        for i = 1:r
          U(:, i) = Ufull(:, inds(n + 1 - i));
          S(i, i) = svals(n + 1 - i);
          V(:, i) = Vfull(:, inds(n + 1 - i));
        end
        
        % Return SVD matrices
        varargout{1} = U;
        varargout{2} = S;
        varargout{3} = V;
    end
  else
    % Sort ascending
    [svals inds] = sort(svals);
    
    % Return full SVD
    if (nargout ~= 3)
        varargout{1} = svals(n:-1:1);
    else
        % Populate full n-dimensional entries
        U = zeros(m);
        S = zeros(m,n);
        V = zeros(n);        
        for i = 1:n
          U(:, i) = Ufull(:, inds(n + 1 - i));
          S(i, i) = svals(n + 1 - i);
          V(:, i) = Vfull(:, inds(n + 1 - i));
        end
        
        %{
        %
        % This code can be used if I ever update mySVD() to compute all m
        % left-singular vectors of a skinny matrix instead of the first n
        % 
        % Poplulate extra (m-n) vectors of U
        nullspaceVectorInds = 1:m;
        nullspaceVectorInds(inds) = [];
        for i = 1:length(nullspaceVectorInds);
          U(:, n+i) = Ufull(:, nullspaceVectorInds(i));
        end
        %}
        
        % Return SVD matrices
        varargout{1} = U;
        varargout{2} = S;
        varargout{3} = V;
    end
  end
end

function dist = SafeDistance(a, b)
  abs_a = abs(a);
  abs_b = abs(b);
  if (abs_a > abs_b)
    dist = abs_a * sqrt(1.0 + abs_b * abs_b / (abs_a * abs_a));
  else
    if (abs_b == 0)
      dist = 0;
    else
      dist = abs_b * sqrt(1.0 + abs_a * abs_a / (abs_b * abs_b));
    end
  end
end
