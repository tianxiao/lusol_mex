classdef lusol < handle
  %lusol  access to LUSOL for computations on a sparse LU factorization.
  %
  % Optional parameters:
  %  options = lusol.luset()
  %  options = lusol.luset('pivot','TRP')
  %
  % Initialization:
  %  lu = lusol(A)
  %
  % Initialization specifying options:
  %  lu = lusol(A,options)
  %
  % Factorize:
  %  [inform nsing depcol] = lu.factorize(A)
  %  [inform nsing depcol] = lu.factorize(A,options)
  %
  % Usage:
  %  y = lu.mulA(x);
  %  x = lu.solveA(b);
  %
  % Update:
  %  inform = lu.repcol(v,j);
  %
  % Table of lusol methods:
  % 
  %   | get options structure   | lusol.luset     |
  %   | main factorize method   | lusol.factorize |
  %   | get factorization stats | lusol.stats     |
  %   | get rank of matrix      | lusol.rank      |
  %   | get inital L factor     | lusol.L0        |
  %   | get U factor            | lusol.U         |
  %   | get row permutation     | lusol.p         |
  %   | get column permutation  | lusol.q         |
  %   | main multiply method    | lusol.mul       |
  %   | compute A*x             | lusol.mulA      |
  %   | compute A'*x            | lusol.mulAt     |
  %   | compute L*x             | lusol.mulL      |
  %   | compute L'*x            | lusol.mulLt     |
  %   | compute U*x             | lusol.mulU      |
  %   | compute U'*x            | lusol.mulUt     |
  %   | main solve method       | lusol.solve     |
  %   | solve A*x = b           | lusol.solveA    |
  %   | solve A'*x = b          | lusol.solveAt   |
  %   | solve L*x = b           | lusol.solveL    |
  %   | solve L'*x = b          | lusol.solveLt   |
  %   | solve U*x = b           | lusol.solveU    |
  %   | solve U'*x = b          | lusol.solveUt   |
  %   | replace a column        | lusol.repcol    |
  %   | replace a row           | lusol.reprow    |
  %   | add a column            | lusol.addcol    |
  %   | add a row               | lusol.addrow    |
  %   | delete a column         | lusol.delcol    |
  %   | delete a row            | lusol.delrow    |
  %   | rank 1 modification     | lusol.r1mod     |
  %
  % See also:
  %  lusol.luset
  %  lusol.factorize
  %  lusol.stats
  %  lusol.L0
  %  lusol.U
  %  lusol.p
  %  lusol.q
  %  lusol.mul
  %  lusol.mulA
  %  lusol.mulAt
  %  lusol.mulL
  %  lusol.mulLt
  %  lusol.mulU
  %  lusol.mulUt
  %  lusol.solve
  %  lusol.solveA
  %  lusol.solveAt
  %  lusol.solveL
  %  lusol.solveLt
  %  lusol.solveU
  %  lusol.solveUt
  %  lusol.repcol
  %  lusol.reprow
  %  lusol.addcol
  %  lusol.addrow
  %  lusol.delcol
  %  lusol.delrow
  %  lusol.r1mod

  properties (Access=private)
    
    % object parameters
    
    minpad = 0; % minimum number of elements to padding arrays
    nzinit = 0; % initial number of non zeros
    minit = 0; % initial number of columns
    ninit = 0; % initial number of rows
    
    % lusol input parameters
    
    maxcol = 0; % max num cols searched for piv element (5)
    pivot = 0; % pivoting method [(0=TPP),1=TRP,2=TCP,3=TSP]
    keepLU = 0; % keep the nonzeros, if 0, permutations are computed (1)
    Ltol1 = 0; % max Lij allowed, default depends on luparm(6)
    Ltol2 = 0; % max Lij allowed during updates
    small = 0; % absolute tolerance for treating reals as zero (eps^0.8)
    Utol1 = 0; % absolute tol for flagging small diags of U (eps^0.67)
    Utol2 = 0; % rel tol for flagging small diags of U (eps^0.67)
    Uspace = 0; % (3.0)
    dens1 = 0; % (0.3)
    dens2 = 0; % (0.5)
    
    % lusol parameter vectors
    
    luparm = 0; % vector of integer parameters (input and output)
    parmlu = 0; % vector of double parameters (input and output)
    
    % scalars
    
    m = 0; % number of rows
    n = 0; % number of columns
    
    mlen = 0; % length of column space vectors
    nlen = 0; % length of row space vectors
    
    nelem = 0; % number of elements in original matrix (may be removed)
    nzmax = 0; % maximum storage allocated for factorization and updates
    diag = 0; % return scalar used by some methods
    vnorm = 0; % return scalar used by some methods
    
    % vectors of lenth nzmax
    
    a = 0; % main storage array
    indc = 0; % row indecies
    indr = 0; % column indecies
    
    % vectors of length mlen
    
    ip = 0; % row permutation
    lenr = 0;
    locr = 0;
    iqloc = 0;
    ipinv = 0;
    v = 0; % storage
    
    % vectors of length nlen
    
    iq = 0; % column permutation
    lenc = 0;
    locc = 0;
    iploc = 0;
    iqinv = 0;
    w = 0; % storage
    
    % other
    
    depcol_lx = 0; % logical index indicating dependent columns after a 
                   % factorize
    
  end
  
  methods (Static)
    function options = luset(varargin)
      %luset  process input to create lusol options structure
      %
      % This method uses Matlab's inputParser class to handle lusol options.
      % It handles no input, input structures, and key-value lists just like
      % many other Matlab programs.
      %
      % To obtain a structure with default settings use:
      %
      %  options = lusol.luset();
      %
      % To create an options structre that will use threshold rook pivoting,
      % use:
      %
      %  options = lusol.luset('pivot','TRP')
      %
      % The best reference for the lusol parameters is currently the comments
      % for the lu1fac subroutine in lusol1.f.
      %
      %
      % |--------+----------+----------------------------------------------------|
      % | param  |  default | description                                        |
      % |--------+----------+----------------------------------------------------|
      %
      % lusol_mex options
      % |--------+----------+----------------------------------------------------|
      % | nzinit |        0 | minimum length for storage arrays                  |
      % | minit  |        0 | minimum length for 'row-space' arrays              |
      % | ninit  |        0 | minimum length for 'column-space' arrays           |
      % | minpad |        5 | minimum padding for row and column space arrays    |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL integer parameters
      % |--------+----------+----------------------------------------------------|
      % | maxcol |        5 | max num cols searched for piv element              |
      % | pivot  |    'TPP' | pivoting method {'TPP','TRP','TCP','TSP'}          |
      % | keepLU |        1 | keep the nonzeros, if 0, permutations are computed |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL real parameters
      % |--------+----------+----------------------------------------------------|
      % | Ltol1  |     10.0 | max Lij allowed, default depends on pivot method   |
      % | Ltol2  |     10.0 | max Lij allowed during updates                     |
      % | small  |  eps^0.8 | absolute tolerance for treating reals as zero      |
      % | Utol1  | eps^0.67 | absolute tol for flagging small diags of U         |
      % | Utol2  | eps^0.67 | rel tol for flagging small diags of U              |
      % | Uspace |      3.0 |                                                    |
      % | dens1  |      0.3 |                                                    |
      % | dens2  |      0.5 |                                                    |
      % |--------+----------+----------------------------------------------------|
      
      
      % get the input parser
      in_parse = inputParser;
      
      % storage parameters
      in_parse.addParamValue('minpad',5,@(x) x>=1);
      in_parse.addParamValue('nzinit',0,@(x) x>=0);
      in_parse.addParamValue('minit',0,@(x) x>=0);
      in_parse.addParamValue('ninit',0,@(x) x>=0);
      
      % lusol integer parameters
      in_parse.addParamValue('maxcol',5,@(x) x>=0);
      in_parse.addParamValue('pivot','TPP',@(x) ismember(x,{'TPP','TRP','TCP','TSP'}));
      in_parse.addParamValue('keepLU',1,@(x) ismember(x,[0 1]));
      
      % lusol real parameters
      in_parse.addParamValue('Ltol1',10.0,@(x) x>=0.0);
      in_parse.addParamValue('Ltol2',10.0,@(x) x>=0.0);
      in_parse.addParamValue('small',eps^0.8,@(x) x>=0.0);
      in_parse.addParamValue('Utol1',eps^0.67,@(x) x>=0.0);
      in_parse.addParamValue('Utol2',eps^0.67,@(x) x>=0.0);
      in_parse.addParamValue('Uspace',3.0,@(x) x>=0.0);
      in_parse.addParamValue('dens1',0.3,@(x) x>=0.0);
      in_parse.addParamValue('dens2',0.5,@(x) x>=0.0);
      
      % parse the input
      in_parse.parse(varargin{:});
      
      % obtain the output
      options = in_parse.Results;
    end
  end
  
  methods (Access=private)
    
    function parse_options(obj,varargin)
      %parse_options  process the options structure to set parameters.

      % use luset to parse user input
      options = obj.luset(varargin{:});
      
      % storage parameters
      obj.minpad = options.minpad;
      obj.nzinit = options.nzinit;
      obj.minit = options.minit;
      obj.ninit = options.ninit;
      
      % lusol integer parameters
      obj.maxcol = options.maxcol;
      obj.keepLU = options.keepLU;
      
      % lusol double parameters
      obj.Ltol1 = options.Ltol1;
      obj.Ltol2 = options.Ltol2;
      obj.small = options.small;
      obj.Utol1 = options.Utol1;
      obj.Utol2 = options.Utol2;
      obj.Uspace = options.Uspace;
      obj.dens1 = options.dens1;
      obj.dens2 = options.dens2;
      
      % set the pivoting strategy
      switch options.pivot
        case 'TPP'
          obj.pivot = 0;
        case 'TRP'
          obj.pivot = 1;
        case 'TCP'
          obj.pivot = 2;
        case 'TSP'
          obj.pivot = 3;
        otherwise
          error('lusol:set_options','Unkown pivot strategy.')
      end
    end
    
    function set_options(obj)
      %set_options allocate and assign parameters to LUSOL arrays
      %
      % LUSOL stores input and output scalar parameters in two vectors:
      %   luparm is an int32 array of length 30
      %   parmlu is a double array of length 30
      %
      % this method sets the LUSOL input parameters in the correct location
      % for the fortran calls.
      %
      
      % allocate parameter vectors
      obj.luparm = zeros(30,1,'int32');
      obj.parmlu = zeros(30,1,'double');
      
      % set parameter values
      obj.luparm(3) = int32(obj.maxcol);
      obj.luparm(6) = int32(obj.pivot);
      obj.luparm(8) = int32(obj.keepLU);
      obj.parmlu(1) = double(obj.Ltol1);
      obj.parmlu(2) = double(obj.Ltol2);
      obj.parmlu(3) = double(obj.small);
      obj.parmlu(4) = double(obj.Utol1);
      obj.parmlu(5) = double(obj.Utol2);
      obj.parmlu(6) = double(obj.Uspace);
      obj.parmlu(7) = double(obj.dens1);
      obj.parmlu(8) = double(obj.dens2);
    end
    
    function allocate(obj,nzmax,mlen,nlen)
      %allocate  allocate LUSOL storage arrays
      %
      % LUSOL operates on many arrays.  This method allocates all of them 
      % to an appropriate size.  Care is taken to not reallocate if it is
      % not needed.
      %
      
      % vectors of length nzmax
      if nzmax > obj.nzmax
        obj.nzmax = int32(nzmax);
        obj.a = zeros(obj.nzmax,1,'double');
        obj.indc = zeros(obj.nzmax,1,'int32');
        obj.indr = zeros(obj.nzmax,1,'int32');
      end
      
      % column space vectors
      % vectors of length mlen
      mpad = mlen + obj.minpad;
      if mlen > obj.mlen
        obj.mlen = int32(mlen);
        obj.ip = zeros(mpad,1,'int32');
        obj.lenr = zeros(mpad,1,'int32');
        obj.locr = zeros(mpad,1,'int32');
        obj.iqloc = zeros(mpad,1,'int32');
        obj.ipinv = zeros(mpad,1,'int32');
        obj.v = zeros(mpad,1,'double');
      end
      
      % vectors of length nlen
      npad = nlen + obj.minpad;
      if nlen > obj.nlen
        obj.nlen = int32(nlen);
        obj.iq = zeros(npad,1,'int32');
        obj.lenc = zeros(npad,1,'int32');
        obj.locc = zeros(npad,1,'int32');
        obj.iploc = zeros(npad,1,'int32');
        obj.iqinv = zeros(npad,1,'int32');
        obj.w = zeros(npad,1,'double');
      end
      
      % other scalars
      obj.diag = double(0.0);
      obj.vnorm = double(0.0);
      
    end
    
    function update_check(obj)
      %update_check  throw and error if this method is called after updates
      %
      % Some methods should not be called after updates to a factorization.
      % This method checks if any updated have occuered and throws and
      % error if this is the case.
      %
      
      nupdat = obj.luparm(15);
      if nupdat > 0
        error('lusol:post_update_call_error','nsing and depcol cannot be called after updates.')
      end
    end
    
  end 
  
  methods
    
    % constructor and main factorize method
    
    function obj = lusol(A,varargin)
      %lusol  constructor for lusol object, factorize A
      %
      % Creates lusol object and factorizes A.  
      %
      % Example:
      %   mylu = lusol(A);
      %
      % See lusol.factorize for more info.
      %
      
      obj.factorize(A,varargin{:});
      
    end
    
    function [inform nsing depcol] = factorize(obj,A,varargin)
      %factorize perform lu factorization on A.
      %
      % This method tells LUSOL to perform an LU factorization on A.
      %
      % Usage (after lu object is initialized):
      %  [inform nsing depcol] = lu.factorize(A)
      %  [inform nsing depcol] = lu.factorize(A,options)
      %
      % Input:
      %  A = matrix to factorize
      %  options = options structure (optional)
      %
      % Output:
      %  inform = status flag
      %  nsing = esimate of the number of singularities
      %  depcol = logical index of dependent columns
      %
      % inform code:
      %  0 if the LU factors were obtained successfully.
      %  1 if U appears to be singular, as judged by lu6chk.
      %  3 if some index pair indc(l), indr(l) lies outside
      %    the matrix dimensions 1:m , 1:n.
      %  4 if some index pair indc(l), indr(l) duplicates
      %    another such pair.
      %  7 if the arrays a, indc, indr were not large enough.
      %    Their length "lena" should be increase to at least
      %    the value "minlen" given in luparm(13).
      %  8 if there was some other fatal error.  (Shouldn't happen!)
      %  9 if no diagonal pivot could be found with TSP or TDP.
      %    The matrix must not be sufficiently definite
      %    or quasi-definite.
      %
      
      % parse optional options
      obj.parse_options(varargin{:});
      obj.set_options();
      
      % get basic information on A
      m = size(A,1);
      n = size(A,2);
      nelem = nnz(A);
      
      if m == 0 || n == 0
        % lusol_mex will cause a segfault if it is given an empty matrix
        % and the pivot method is TRP.  It is best to just not use LUSOL on
        % an empty matrix.
        error('lusol:factorize','LUSOL cannot factorize an empty matrix.')
      end
      
      % set storage sizes
      nzmax = max([2*nelem 10*m 10*n 10000 obj.nzinit]);
      mlen = max(m,obj.minit);
      nlen = max(n,obj.ninit);
      
      % nzmax, mlen, and nlen are all now large enough to load the matrix
      % into the appropriate arrays.  However, nzmax may be too small to
      % factorize.  In this case LUSOL should return an error.
      
      % possibly (re)allocated storage arrays
      obj.allocate(nzmax,mlen,nlen);
      
      % set size and number of elements appropriatly
      obj.nelem = int32(nelem);
      obj.m = int32(m);
      obj.n = int32(n);
      
      % extract data from A for use in LUSOL
      [i j s] = find(A);
      obj.indc(1:obj.nelem) = int32(i);
      obj.indr(1:obj.nelem) = int32(j);
      obj.a(1:obj.nelem) = double(s);
      
      % run lusol
      ret_inform = int32(0);
      lusol_mex('lu1fac', ...
        obj.m, ...
        obj.n, ...
        obj.nelem, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        obj.iploc, ...
        obj.iqloc, ...
        obj.ipinv, ...
        obj.iqinv, ...
        obj.w, ...
        ret_inform);
      
      % error checking
      switch ret_inform
        case 0
          % ok, LU factors obtained
        case 1
          % ok, LU factors obtained, rank deficient
        case 3
          % error, some index pair indc(l), indr(l) lies outside
          % the matrix dimensions 1:m , 1:n.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize','LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 4
          % error, some index pair indc(l), indr(l) duplicates
          % another such pair.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize','LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 7
          % error, not enough storage.  User needs to increase nzinit
          % parameter
          err = MException('lusol:factorize','LUSOL needs more storage.  Increase the nzinit parameter.  inform = %d',ret_inform);
          throw(err);
        case 8
          % some other fatal error
          err = MException('lusol:factorize','LUSOL other fatal error. inform = %d',ret_inform);
          throw(err);
        case 9
          % error, no diagonal pivot could be found with TSP or TDP.
          % The matrix must not be sufficiently definite or quasi-definite
          err = MException('lusol:factorize','LUSOL no diagonal pivot could be found with TSP or TDP. inform = %d',ret_inform);
          throw(err);
      end
      
      obj.depcol_lx = (obj.w(1:obj.n) <= 0.0);
      
      % user requests inform flag
      if nargout > 0
        inform = obj.inform();
      end
      
      % user requests number of singularities
      if nargout > 1
        nsing = obj.nsing();
      end
      
      % user requests dependent column indicator
      if nargout > 2
        depcol = obj.depcol_lx;
      end
      
    end
    
    % methods to collect information about matrix and factorization
    
    function [m n] = size(obj)
      %size  get size of factorized matrix
      m = double(obj.m);
      n = double(obj.n);
    end
    
    function s = stats(obj)
      %stats  return LUSOL stats structure
      %
      % This method builds a Matlab struct containing LUSOL output
      % parameters.
      %
      % LUSOL output parameters:
      %
      % inform   Return code from last call to any LU routine.
      % nsing    No. of singularities marked in the
      %          output array w(*).
      % jsing    Column index of last singularity.
      % minlen   Minimum recommended value for  lena.
      % maxlen   ?
      % nupdat   No. of updates performed by the lu8 routines.
      % nrank    No. of nonempty rows of U.
      % ndens1   No. of columns remaining when the density of
      %          the matrix being factorized reached dens1.
      % ndens2   No. of columns remaining when the density of
      %          the matrix being factorized reached dens2.
      % jumin    The column index associated with DUmin.
      % numL0    No. of columns in initial  L.
      % lenL0    Size of initial  L  (no. of nonzeros).
      % lenU0    Size of initial  U.
      % lenL     Size of current  L.
      % lenU     Size of current  U.
      % lrow     Length of row file.
      % ncp      No. of compressions of LU data structures.
      % mersum   lu1fac: sum of Markowitz merit counts.
      % nUtri    lu1fac: triangular rows in U.
      % nLtri    lu1fac: triangular rows in L.
      % Amax     Maximum element in  A.
      % Lmax     Maximum multiplier in current  L.
      % Umax     Maximum element in current  U.
      % DUmax    Maximum diagonal in  U.
      % DUmin    Minimum diagonal in  U.
      % Akmax    Maximum element generated at any stage
      %          during TCP factorization.
      % growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax.
      % resid    lu6sol: residual after solve with U or U'.
      
      s.inform = double(obj.luparm(10));
      s.nsing = double(obj.luparm(11));
      s.jsing = double(obj.luparm(12));
      s.minlen = double(obj.luparm(13));
      s.maxlen = double(obj.luparm(14));
      s.nupdat = double(obj.luparm(15));
      s.nrank = double(obj.luparm(16));
      s.ndens1 = double(obj.luparm(17));
      s.ndens2 = double(obj.luparm(18));
      s.jumin = double(obj.luparm(19));
      s.numL0 = double(obj.luparm(20));
      s.lenL0 = double(obj.luparm(21));
      s.lenU0 = double(obj.luparm(22));
      s.lenL = double(obj.luparm(23));
      s.lenU = double(obj.luparm(24));
      s.lrow = double(obj.luparm(25));
      s.ncp = double(obj.luparm(26));
      s.mersum = double(obj.luparm(27));
      s.nUtri = double(obj.luparm(28));
      s.nLtri = double(obj.luparm(29));
      s.Amax = double(obj.parmlu(10));
      s.Lmax = double(obj.parmlu(11));
      s.Umax = double(obj.parmlu(12));
      s.DUmax = double(obj.parmlu(13));
      s.DUmin = double(obj.parmlu(14));
      s.Akmax = double(obj.parmlu(15));
      s.growth = double(obj.parmlu(16));
      s.resid = double(obj.parmlu(20));
      
    end
    
    function info = inform(obj)
      %inform  return code from last call to LUSOL routines
      info = double(obj.luparm(10));
    end
    
    function k = nsing(obj)
      %nsing  number of singularities marked in depcol
      %
      % This method may not be called after updates.  The nsing parameter
      % is only computed after a full factorize.
      %
      
      % this method only works if no updates have occured.
      obj.update_check();
      
      k = double(obj.luparm(11));
    end
    
    function k = rank(obj)
      %rank  the rank of the matrix determined by the number of independent columns
      %
      % This method uses the LUSOL parameter nsing.  This is only computed
      % after a full factorize.  Thus this method should not be used to
      % determine the rank of a matrix after updates.  In that case look at
      % the flags that are returned by the update methods.
      %
      % Rank determination with LUSOL is more reliable under threshold rook
      % pivoting.  Example options:
      %
      %  options = lusol.luset('pivot','TRP','Ltol1',10)
      %
      
      k = double(obj.n) - double(obj.luparm(11));
    end
    
    function d = depcol(obj)
      %depcol  logical vector indicating dependent columns
      %
      % This method may not be called after updates.  The method looks at
      % data that only relevant after a factorize.
      %
      
      % this method only works if no updated have occured.
      obj.update_check();
      
      d = obj.depcol_lx;
    end
    
    function ip = p(obj)
      %p  return row permutation vector
      ip = double(obj.ip(1:obj.m));
    end
    
    function iq = q(obj)
      %q  return column permutation vector
      iq = double(obj.iq(1:obj.n));
    end
    
    % methods to get the matrix factors
    
    function [U p q] = U(obj,pm_opt)
      %U  get the upper triangular factor U
      %
      % Extract the U factor from LUSOL data and return as a Matlab sparse
      % matrix.
      %
      % Set up:
      %   mylu = lusol(A);
      %
      % Usage:
      %   U1 = mylu.U();
      %   [U2 p q] = mylu.U();
      %   [U2 P Q] = mylu.U('matrix');
      %
      % The first call returns U as a permuted triangle.  The second
      % returns U as an upper triangular matrix with permutation vectors p
      % and q.  The third call returns sparse permutation matrices P and Q.
      % The result of the three calls would produce:
      %   U1(p,1) == P*U1*Q == U2 == upper triangular
      %
      
      %
      % After a factorize (call to lu1fac) LUSOL stores U by rows at the
      % start of arrays a and indr.  lenr(1:m) stores the number of entries
      % in each row in original order.  locr(1:m) points to the beginning
      % of rows and is stored in original order.
      %
      % Special care must be taken when A is rank deficient.  LUSOL
      % actually stores lenU-nsing entries.  I suppose the extra nsing
      % contained in lenU could be for the zeros on the diagonal.  However,
      % LUSOL seems to handle these implicitly.
      %
      
      % permutation flag, set true if user desires upper triangular U and
      % permutation vectors
      permflg = false;
      
      % matrix flag, set true if user desires sparse permutation matrices
      % instead of vectors
      matrflg = false;
      
      if nargout == 3
        permflg = true;
      end
      
      if nargin >= 2 && strcmp(pm_opt,'matrix')
        matrflg = true;
      end
      
      [m n] = obj.size();
      p = obj.p();
      q = obj.q();
      s = obj.stats();
      
      % initialize arrays for U triplets
      ui = zeros(s.lenU-s.nsing,1,'double');
      uj = zeros(s.lenU-s.nsing,1,'double');
      ua = zeros(s.lenU-s.nsing,1,'double');
      
      % array position pointers
      k1 = 1;
      k2 = 1;
      
      for i = 1:s.nrank
        % get row index
        piv = p(i);
        
        % get length of row
        len = double(obj.lenr(piv));
        
        % get location of row
        loc = double(obj.locr(piv));
        
        k2 = k1+len-1;
        
        % load data into triplet arrays
        ui(k1:k2) = piv*ones(len,1);
        uj(k1:k2) = double(obj.indr(loc:loc+len-1));
        ua(k1:k2) = double(obj.a(loc:loc+len-1));
        
        k1 = k1+len;
        
      end
      
      % generate sparse matrix
      U = sparse(ui,uj,ua,m,n);
      
      if permflg
        % produce and return upper triangular U
        U = U(p,q);
      end
      
      if matrflg
        % construct and return sparse permutation matrices
        p = sparse(1:n,p,1,n,n);
        q = sparse(q,1:m,1,m,m);
      end
      
    end
    
    function [L0 p] = L0(obj,pm_opt)
      %L0  get the initial lower triangular factor L0
      %
      % Extracts the initial lower triangular factor from the LUSOL data
      % structure.  LUSOL stores updates to L in product form, thus updates
      % are not included in L0.
      %
      % Set up:
      %   mylu.factorize(A);
      %
      % Usage:
      %   L1 = mylu.get_L0();
      %   [L2 p] = mylu.get_L0();
      %   [L2 P] = mylu.get_L0('matrix');
      %
      % The first call returns L1 as a permuted triangle.  The second and
      % third call will return L2 as a lower triangular matrix with
      % permutation vector p or matrix P.  The result will give:
      %   L1(p,p) == P*L1*P' == L2 == lower triangular
      %
      
      %
      % After a factorize (call to lu1fac) LUSOL stores non-trivial columns
      % of L at the end of a, indc, and indr.  lenc(1:numL0) stores the
      % number of entries in each column, not including the 1 on the
      % diagonal.  The negatives of the elements of L are stored in a.
      %
      % indc gives the row indices for non-zero elements
      % indr gives the column indices
      %
      % It turns out that off diagonal elements of L are stored in triplet
      % form at the end of a, indc, and indr.  It remains to add the ones
      % on the diagonal.
      %
      
      % permutation flag, set true if user desires upper triangular U and
      % permutation vectors
      permflg = false;
      
      % matrix flag, set true if user desires sparse permutation matrices
      % instead of vectors
      matrflg = false;
      
      if nargout == 2
        permflg = true;
      end
      
      if nargin >= 2 && strcmp(pm_opt,'matrix')
        matrflg = true;
      end
      
      % obtain information from object
      [m n] = obj.size();
      p = obj.p();
      s = obj.stats();
      lena = double(obj.nzmax);
      
      % allocate arrays for triplet form
      li = zeros(s.lenL0+m,1);
      lj = zeros(s.lenL0+m,1);
      la = zeros(s.lenL0+m,1);
      
      % read the triplet form from LUSOL data
      li(1:s.lenL0) = double(obj.indc(lena-s.lenL0+1:lena));
      lj(1:s.lenL0) = double(obj.indr(lena-s.lenL0+1:lena));
      la(1:s.lenL0) = -double(obj.a(lena-s.lenL0+1:lena));
      
      % add 1's along the diagonal
      li(s.lenL0+1:end) = (1:m)';
      lj(s.lenL0+1:end) = (1:m)';
      la(s.lenL0+1:end) = ones(m,1);
      
      % create matlab sparse matrix
      L0 = sparse(li,lj,la);
      
      if permflg
        % produce and return lower triangular L0
        L0 = L0(p,p);
      end
      
      if matrflg
        % construct and return sparse permutation matrix
        p = sparse(1:n,p,1);
      end
      
    end
    
    % solve methods
    
    function [x inform resid] = solve(obj,b,mode)
      %solve  call lu6sol to perform various solves with L and U factors.
      %
      % Usage:
      %  x = lu.solve(b,mode)
      %
      % Input:
      %  b = right hand side vector
      %  mode = solution mode (see table below)
      %
      % Output:
      %  x = solution vector
      %  inform = status flag
      %  resid = 1-norm of residual
      %
      % Modes:
      %  1    x  solves   L x = b
      %  2    x  solves   L'x = b
      %  3    x  solves   U x = b
      %  4    x  solves   U'x = b
      %  5    x  solves   A x = b (default)
      %  6    x  solves   A'x = b
      %
      % inform flags:
      %  0 = successful solve
      %  1 = if U is singular, and residual is non-zero
      %
      
      if nargin < 3
        mode = 5;
      end
      
      mode = int32(mode);
      
      % force matlab to copy
      % I think I can remove this line because I force a copy in the switch
      % b(1) = b(1);
      
      % make sure b is a vector
      if ~isvector(b)
        error('lusol:solve','b must be a vector.')
      end
      
      % orient b vector
      b = double(full(b(:)));
      lenb = length(b);
      
      switch mode
        case 1
          if lenb ~= obj.m, error('lusol:solve','b has incorrect size.'); end
          obj.v(1:obj.m) = b;
        case 2
          if lenb ~= obj.m, error('lusol:solve','b has incorrect size.'); end
          obj.v(1:obj.m) = b;
        case 3
          if lenb ~= obj.m, error('lusol:solve','b has incorrect size.'); end
          obj.v(1:obj.m) = b;
        case 4
          if lenb ~= obj.n, error('lusol:solve','b has incorrect size.'); end
          obj.w(1:obj.n) = b;
        case 5
          if lenb ~= obj.m, error('lusol:solve','b has incorrect size.'); end
          obj.v(1:obj.m) = b;
        case 6
          if lenb ~= obj.n, error('lusol:solve','b has incorrect size.'); end
          obj.w(1:obj.n) = b;
        otherwise
          error('lusol:solve','unrecognized mode.')
      end
      
      ret_inform = int32(0);
      lusol_mex('lu6sol', ...
        mode, ...
        obj.m, ...
        obj.n, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform);
      
      switch mode
        case 1
          x = obj.v(1:obj.m);
        case 2
          x = obj.v(1:obj.m);
        case 3
          x = obj.w(1:obj.n);
        case 4
          x = obj.v(1:obj.m);
        case 5
          x = obj.w(1:obj.n);
        case 6
          x = obj.v(1:obj.m);
      end
      
      inform = obj.luparm(10);
      resid = obj.parmlu(20);
    end
    function [x inform resid] = solveA(obj,b)
      %solveA  solve A x = b.
      %
      % see also lusol.solve
      [x inform resid] = obj.solve(b,5);
    end
    function [x inform resid] = solveAt(obj,b)
      %solveAt  solve A'x = b.
      %
      % see also lusol.solve
      [x inform resid] = obj.solve(b,6);
    end
    function [x inform] = solveL(obj,b)
      %solveL  solve L x = b.
      %
      % see also lusol.solve
      [x inform] = obj.solve(b,1);
    end
    function [x inform] = solveLt(obj,b)
      %solveLt  solve L'x = b.
      %
      % see also lusol.solve
      [x inform] = obj.solve(b,2);
    end
    function [x inform resid] = solveU(obj,b)
      %solveU  solve U x = b.
      %
      % see also lusol.solve
      [x inform resid] = obj.solve(b,3);
    end
    function [x inform resid] = solveUt(obj,b)
      %solveUt  solve U'x = b.
      %
      % see also lusol.solve
      [x inform resid] = obj.solve(b,4);
    end
    
    % multiply methods
    
    function y = mul(obj,x,mode)
      %mul  call LUSOL to perform various multiplies with L and U factors.
      %
      % Usage:
      %  y = lu.mul(x,mode)
      %
      % mode
      %  1    y = L x
      %  2    y = L'x
      %  3    y = U x
      %  4    y = U'x
      %  5    y = A x (default)
      %  6    y = A'x
      %
      % Warning: it seems like mulA works, but mulAt does not.
      %
      if nargin < 3
        mode = 5;
      end
      
      mode = int32(mode);
      
      % check if x is a vector
      if ~isvector(x)
        error('lusol:mul','x must be a vector.');
      end
      
      % force matlab to copy
      % not needed with padding method
      % x(1) = x(1);
      
      % orient x vector
      x = double(full(x(:)));
      lenx = length(x);
      
      switch mode
        case 1
          if lenx ~= obj.m, error('lusol:mul','x has incorrect size.'); end
          obj.v(1:obj.m) = x;
        case 2
          if lenx ~= obj.m, error('lusol:mul','x has incorrect size.'); end
          obj.v(1:obj.m) = x;
        case 3
          if lenx ~= obj.n, error('lusol:mul','x has incorrect size.'); end
          obj.w(1:obj.n) = x;
        case 4
          if lenx ~= obj.m, error('lusol:mul','x has incorrect size.'); end
          obj.v(1:obj.m) = x;
        case 5
          if lenx ~= obj.n, error('lusol:mul','x has incorrect size.'); end
          obj.w(1:obj.n) = x;
        case 6
          if lenx ~= obj.m, error('lusol:mul','x has incorrect size.'); end
          obj.v(1:obj.m) = x;
        otherwise
          error('lusol:mul','unrecognized mode.')
      end
      
      lusol_mex('lu6mul', ...
        mode, ...
        obj.m, ...
        obj.n, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr);
      
      switch mode
        case 1
          y = obj.v(1:obj.m);
        case 2
          y = obj.v(1:obj.m);
        case 3
          y = obj.v(1:obj.m);
        case 4
          y = obj.w(1:obj.n);
        case 5
          y = obj.v(1:obj.m);
        case 6
          y = obj.w(1:obj.n);
      end
      
    end
    function y = mulA(obj,x)
      %mulA  compute y = A x.
      %
      % see also lusol.mul
      y = obj.mul(x,5);
    end
    function y = mulAt(obj,x)
      %mulAt  compute y = A'x.
      %
      % Warning: this does not seem to work at the moment.
      %
      % see also lusol.mul
      y = obj.mul(x,6);
    end
    function y = mulL(obj,x)
      %mulL  compute y = L x.
      %
      % see also lusol.mul
      y = obj.mul(x,1);
    end
    function y = mulLt(obj,x)
      %mulLt  compute y = L'x.
      %
      % see also lusol.mul
      y = obj.mul(x,2);
    end
    function y = mulU(obj,x)
      %mulU  compute y = U x.
      %
      % see also lusol.mul
      y = obj.mul(x,3);
    end
    function y = mulUt(obj,x)
      %mulUt  compute y = U'x.
      %
      % see also lusol.mul
      y = obj.mul(x,4);
    end
    
    % update methods
    
    function [inform diag vnorm] = repcol(obj,v,j)
      %repcol  update LU factorization to replace a column
      %
      % Usage:
      %  [inform diag vnorm] = lu.repcol(v,j)
      %
      % Inputs:
      %  v = new column
      %  j = column to replace
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %  vnorm = ?
      %
      % On exit:
      %  inform = -1  if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  2  if the update seemed to be unstable
      %               (diag much bigger than vnorm).
      %  inform =  7  if the update was not completed (lack of storage).
      %  inform =  8  if j is not between 1 and n.
      %
      
      if ~isvector(v) || length(v) ~= obj.m
        error('lusol:repcol','v must be a vector of length m.')
      end
      
      % force copy and orient
      % not needed with pad method
      % v(1) = v(1);
      % v = v(:);
      
      % densify and copy
      obj.v(1:obj.m) = double(full(v(:)));
      
      mode1 = int32(1);
      mode2 = int32(1);
      j = int32(j);
      
      ret_inform = int32(0);
      lusol_mex('lu8rpc', ...
        mode1, ...
        mode2, ...
        obj.m, ...
        obj.n, ...
        j, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform, ...
        obj.diag, ...
        obj.vnorm);
      
      inform = obj.inform();
      diag = obj.diag;
      vnorm = obj.vnorm;
    end
    
    function inform = reprow(obj,w,i)
      %reprow  update LU factorization to replace a row
      %
      % Usage:
      %  inform = lu.reprow(w,i)
      %
      % Inputs:
      %  w = new row
      %  i = row to replace
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform = -1  if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %  inform =  8  if i is not between 1 and m.
      %
      
      if ~isvector(w) || length(w) ~= obj.n
        error('lusol:reprow','w must be a vector of length n.')
      end
      
      % force copy and orient
      % not needed with padding method
      % w(1) = w(1);
      % w = w(:);
      
      % densify and copy
      wnew = zeros(obj.nlen,1,'double');
      wnew(1:obj.nlen) = double(full(w(:)));
      
      mode1 = int32(1);
      mode2 = int32(1);
      i = int32(i);
      
      ret_inform = int32(0);
      lusol_mex('lu8rpr', ...
        mode1, ...
        mode2, ...
        obj.m, ...
        obj.n, ...
        i, ...
        obj.v, ...
        obj.w, ...
        wnew, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform);
      
      inform = obj.luparm(10);
      
    end
    
    function [inform diag vnorm] = addcol(obj,v)
      %addcol  update LU factorization to add column to end of A
      %
      % Usage:
      %  [inform diag vnorm] = lu.addcol(v)
      %
      % Input:
      %  v = vector of length m to added to the end of A
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %  vnorm = ?
      %
      % On exit:
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % Note that this will change n and the size of some workspace vectors.
      %
      
      % make sure n < nlen-minpad
      if obj.n == obj.nlen-obj.minpad
        error('lusol:addcol','Not enough space to add a column.  Use larger nlen.')
      end
      
      % check input
      if ~isvector(v) || length(v) ~= obj.m
        error('lusol:addcol','v must be a vector of length m.')
      end
      
      % force copy and orient
      % v(1) = v(1);
      % v = v(:);
      
      % densify
      obj.v(1:obj.m) = double(full(v(:)));
      
      % set mode
      mode = int32(1);
      
      % increment n
      obj.n = int32(obj.n + 1);
      
      % set values in nlen vectors
      obj.w(obj.n) = 0.0;
      obj.iq(obj.n) = obj.n;
      obj.lenc(obj.n) = int32(0);
      obj.locc(obj.n) = int32(0);
      
      ret_inform = int32(0);
      lusol_mex('lu8adc', ...
        mode, ...
        obj.m, ...
        obj.n, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform, ...
        obj.diag, ...
        obj.vnorm);
      
      inform = obj.luparm(10);
      diag = obj.diag;
      vnorm = obj.vnorm;
    end
    
    function [inform diag] = addrow(obj,w)
      %addrow  update LU factorization to add row to end of A
      %
      % Usage:
      %  [inform diag] = lu.addrow(w)
      %
      % Input:
      %  w = vector of length n to added to the bottom of A
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %
      % On exit:
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This will change m and the size of some workspace vectors.
      %
      
      % make sure m < mlen-minpad
      if obj.m == obj.mlen-obj.minpad
        error('lusol:addrow','Not enough space to add a row.  Use larger mlen.')
      end
      
      % check input
      if ~isvector(w) || length(w) ~= obj.n
        error('lusol:addrow','w must be a vector of length n.')
      end
      
      % force copy and orient
      % w(1) = w(1);
      % w = w(:);
      
      % densify
      obj.w(1:obj.n) = double(full(w(:)));
      
      % increment m
      obj.m = int32(obj.m + 1);
      
      % set values in mlen vectors
      obj.v(obj.m) = 0.0;
      obj.ip(obj.m) = obj.m;
      obj.lenr(obj.m) = int32(0);
      obj.locr(obj.m) = int32(0);
      
      ret_inform = int32(0);
      lusol_mex('lu8adr', ...
        obj.m, ...
        obj.n, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform, ...
        obj.diag);
      
      inform = obj.luparm(10);
      diag = obj.diag;
      
    end
    
    function inform = delcol(obj,j)
      %delcol  update LU factorization to delete column j from A
      %
      % Usage:
      %  inform = lu.delcol(j)
      %
      % Input:
      %  j = column to delete from A
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This will decrement n by 1 and change relevant storage vectors.
      %
      
      % check to make sure there is more than 1 column
      if obj.n == 1
        error('lusol:delcol','matrix has only one column, cannot delete.')
      end
      
      % check if j is within bounds
      if j < 1 || j > obj.n
        error('lusol:delcol','j is out of bounds.');
      end
      
      j = int32(j);
      
      ret_inform = int32(0);
      lusol_mex('lu8dlc', ...
        obj.m, ...
        obj.n, ...
        j, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform);
      
      % decrement n
      obj.n = int32(obj.n-1);
      
      inform = obj.luparm(10);
    end
    
    function inform = delrow(obj,i)
      %delrow  update LU factorization to delete row i from A
      %
      % Usage:
      %  inform = lu.delrow(i)
      %
      % Input:
      %  i = row to delete from A
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This does not change the size of A and leaves zeros in the last row.
      %
      
      % check to make sure there is more than 1 row
      if obj.m == 1
        error('lusol:delrow','matrix has only one row, cannot delete.')
      end
      
      % check input
      if i < 1 || i > obj.m
        error(lusol:delrow,'i is out of bounds.');
      end
      
      i = int32(i);
      
      % set mode
      % If  mode = 1,  the old row is assumed to be unknown.  It will be
      %                computed from the LU factors of  A.
      % If  mode = 2,  w  must contain the old row.
      mode = int32(1);
      
      ret_inform = int32(0);
      lusol_mex('lu8dlr', ...
        mode, ...
        obj.m, ...
        obj.n, ...
        i, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform);
      
      % from LUSOL code comments.  This method does not change the size of
      % A.  It introduces a zero row, then permutes it to the bottom.  This
      % will have to be tested.  For now, I am not going to change the
      % value of m.
      
      inform = obj.luparm(10);
      
    end
    
    function inform = r1mod(obj,l,u,beta)
      %r1mod  update LU factorization to perform rank-1 update A+beta*l*u'
      %
      % Usage:
      %  inform = lu.r1mod(l,u,beta)
      %
      % Input:
      %  l = vector of length m
      %  u = vector of length n
      %  beta = scalar value
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      
      % handle variable input
      if nargin < 4 || isempty(beta)
        beta = 1.0;
      end
      
      % check input
      if ~isvector(l) || length(l) ~= obj.m
        error('lusol:r1mod','l must be a vector of length m.')
      end
      if ~isvector(u) || length(u) ~= obj.n
        error('lusol:r1mod','u must be a vector of length n.')
      end
      if ~isscalar(beta)
        error('lusol:r1mod','beta must be a scalar.')
      end
      
      % orient, cast, and copy
      obj.v(1:obj.mlen) = double(full(l(:)));
      obj.w(1:obj.nlen) = double(full(u(:)));
      beta = double(beta);
      
      % set mode
      mode = int32(1);
      
      ret_inform = int32(0);
      lusol_mex('lu8mod', ...
        mode, ...
        obj.m, ...
        obj.n, ...
        beta, ...
        obj.v, ...
        obj.w, ...
        obj.nzmax, ...
        obj.luparm, ...
        obj.parmlu, ...
        obj.a, ...
        obj.indc, ...
        obj.indr, ...
        obj.ip, ...
        obj.iq, ...
        obj.lenc, ...
        obj.lenr, ...
        obj.locc, ...
        obj.locr, ...
        ret_inform);
      
      inform = obj.luparm(10);
      
    end
    
  end 
  
end % classdef