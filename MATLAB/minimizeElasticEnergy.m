function newVertex = minimizeElasticEnergy(face, vertex, tarLength, ...
    varargin )
%MINIMIZEELASTICENERGY this function minimizes the elastic energy of a
%non-Euclidean shell with a given initial configuration and a specified
%intrinsic geometry. This is the primary MATLAB interface for the
%(N)on-(E)uclidean (S)hell (S)imulator (NES).
%
%   INPUT PARAMETERS:
%
%       - face:         #Fx3 face connectivity list
%
%       - vertex:       #Vx3 3D vertex coordinate list
%
%       - tarLength:    #Ex1 list of target edge lengths
%
%   OPTIONAL INPUT PARAMETERS: (Name, Value)-Pairs
%
%       - ('TargetAngles', tarTheta = zeros(size(tarLength))):
%           #Ex1 list of target edge bending angles (i.e. the angle
%           between the normal vectors of adjacent faces). Undefined for
%           boundary edges (but should be set to zero). Default behavior is
%           that of a non-Euclidean plate
%
%       - ('Thickness', h = 0.01):
%           The thickness of the elastic sheet in the same units as the
%           edge lengths
%
%       - ('Poisson', nu = 0.5):
%           The Poisson ratio of the elastic sheet.
%
%       - ('Alpha', alpha = 0):
%           The scalar coefficient of the fixed point energy. Set to zero
%           if no vertices are pinned to target locations
%
%       - ('TargetVertices', target_ID = []):
%           The vertex IDs of any vertices that are pinned to target
%           locations
%
%       - ('TargetLocations', targetLocations = []):
%           The 3D target locations of the fixed vertices. Does not have to
%           be equal to the initial locations of the vertices
%
%       - ('FixBoundary'):
%           If this option is flagged, all boundary vertices of the input
%           mesh will be fixed to their initial locations.
%           (Default alpha = 1)
%
%       - ('Beta', beta = 0):
%           The scalar coefficient of the fixed volume energy. Set to zero
%           if no target volume is desired.
%
%       - ('TargetVolume', targetVolume = []):
%           The target volume enclosed by the elastic body
%
%       - ('FixVolume'):
%           If this option is flagged, the volume is fixed to the initial
%           volume of the elastic body (Default beta = 1)
%
%       - ('PhantomFaces', phantomFaces = []):
%           An alternative (#F'x3) face connectivity list that fills any
%           holes in the surface, but leaves the vertex number and
%           locations unchanged. This allows the user to set a target
%           volume for surfaces with boundaries.
%
%       - ('Mu', mu = 0):
%           The inverse scalar coefficient of the growth restriction
%           energy. Set to zero if growth is not restricted in any way
%
%       - ('GrowthRestrictionField', restrictVector = []):
%           A (#Fx3) set of tangent unit vectors to the surface defined on
%           mesh faces. If supplied, edge vectors in the face will be
%           constrained to not exceed a maximum projected length along the
%           restriction vector.
%
%       - ('RestrictedLengths', restrictLengths = []):
%           A (#Fx3) set of maximum projected edge lengths along the growth
%           restriction vector for each edge in a given face
%
%       - ('MHistorySize', m = 20):
%           An L-BFGS parameter. The number of corrections used to
%           approximate the inverse Hessian matrix. The L-BFGS routine
%           stores the computation results of the previous m iterations to
%           approximate the inverse Hessian matrix of the current
%           iteration. This parameter controls the size of the limited
%           memories (corrections). Values less than 3 are not recommended.
%           Large values can result in excessive computing time
%
%       - ('Epsilon', epsilon = 1e-5):
%           An L-BFGS parameter. Absolute tolerance for the gradient norm
%           convergence test. A minimization terminates when the L2 norm of
%           the gradient ||g|| < max(epsilon, epsilon_rel * ||x||), where x
%           is the vector of independent variables
%
%       - ('EpsilonRel', epsilon_rel = 1e-7):
%           An L-BFGS parameter. Relative tolerance for the gradient norm
%           convergence test. A minimization terminates when the L2 norm of
%           the gradient ||g|| < max(epsilon, epsilon_rel * ||x||), where x
%           is the vector of independent variables
%
%       - ('Past', past = 0):
%           An L-BFGS parameter. Distance for delta-based convergence test.
%           This parameter determines the number of iterations d over which
%           to compute the rate of decrease of the objective function
%           f_{k-d}(x)-f_k(x), where k is the current iteration. If the
%           value of this parameters is zero, the delta-based convergence
%           criterion is ignored
%
%       - ('Delta', delta = 0):
%           An L-BFGS parameter. Threshold for delta-based convergence
%           test. A minimization terminates when
%           |f_{k-d}(x)-f_k(x)| < delta * max(1, |f_k(x)|, |f_{k-d}(x)|)
%           where f_k(x) is the current function value and f_{k-d}(x) is
%           the function value d iterations ago
%
%       - ('MaxIterations', max_iterations = 10000):
%           An L-BFGS parameter. The maximum number of allowed iterations.
%           A mimimization terminates when the iteration count exceeds this
%           number. Setting this parameter to zero continues a minimization
%           procedure until convergence or an error
%
%       - ('MaxLinesearch', max_linesearch = 20):
%           An L-BFGS parameter. The maximum number of trials for the line
%           search. The parameter controls the number of function and
%           gradient evaluations per iteration of the line search routine
%
%       - ('MinStep', min_step = 1e-20):
%           An L-BFGS parameter. The minimum step length allowed in the
%           line search. Usually this value does not need to be modified
%
%       - ('MaxStep', max_step = 1e20):
%           An L-BFGS parameter. The maximum step length allowd in the line
%           search. Usually this value does not need to be modified
%
%       - ('FTol', ftol = 1e-4):
%           An L-BFGS parameter. A parameter to control the accuracy of the
%           line search routin. This parameter should be greater than zero
%           and smaller than 0.5. Usually this value does not need to be
%           modified
%
%       - ('WolfeCoeff', wolfe = 0.9):
%           An L-BFGS parameter. The coefficient for the Wolfe condition.
%           This parameter is only valid when the line search algorithm is
%           used with the Wolfe (or strong Wolfe) termination conditions.
%           This parameter should be greater than the 'ftol' parameter and
%           smaller than 1.0. Usually this value does not need to be
%           modified
%
%       - ('IterDisplay', iterDisp = false);
%           An L-BFGS parameter. If true, the mimization will display
%           iterative updates and a detailed description of the
%           optimization termination condition
%
%       - ('Linesearch', 'StrongWolfe'):
%           An L-BFGS parameter. The line search termination condition
%
%           (1) 'Armijo'. The line search method finds a step length that
%           satisfies the sufficient decrease (Armijo) condition,
%           f(x+a*d) <= f(x) + ftol * a * g(x)^T * d
%           where x is the current point, d is the current search
%           direction, a is the step length, and f(x) and g(x) are teh
%           function and gradient values respectively
%
%           (2) 'Wolfe'. The line search method finds a step length that
%           satisfies BOTH the Armijo condition and the curvature condition
%           g(x+a*d)^T * d >= wolfe * g(x)^T * d
%
%           (3) 'StrongWolfe'. The line search method finds a step length
%           that satisfies BOTH the Armijo condtion and the strong
%           curvature condition, |g(x+a*d)^T * d| >= wolfe * |g(x)^T * d|
%          
%       - ('LinesearchMethod', 'Backtracking'):
%           An L-BFGS parameter. The line search algorithm
%
%           (1) 'Backtracking'. A line search method where a large initial
%           guess is tried and then backtracked until a suitable candidate
%           is found. Compatible with all line search termination
%           conditions NOTE: This is the only linesearch method that is
%           currently safe against Inf/NaN values in the objective function
%           and gradient
%
%           (2) 'Bracketing'. Similar to the backtracking line search
%           method except that it actively maintains an upper and lower
%           bound of the current search range. Compatible with all line
%           search termination conditions
%
%           (3) 'NocedalWright'. A line search algorithm that guarantees to
%           find a step size that satisfies the strong Wolfe conditions.
%           Implementation is based on
%           "Numerical Optimzation" 2nd Edition,
%           Jorge Nocedal Stephen J. Wright
%           Chapter 3. Line Search Methods, pg 60
%
%   by Dillon Cislo
%   Copyright (C) 2019-2021 Dillon Cislo <dilloncislo@gmail.com>
%   Under GPL license

%--------------------------------------------------------------------------
% Set Default Parameter Values
%--------------------------------------------------------------------------

% Material parameters -----------------------------------------------------
h = 0.01;
nu = 0.5;

% Target vertex correspondence parameters ---------------------------------
alpha = 0;
target_ID = [];
targetLocations = [];

fixBoundary = false;
anyFixed = false;

% Target volume correspondence parameters ---------------------------------
beta = 0;
targetVolume = [];
phantomFaces = [];

fixVolume = false;
usePhantom = false;

% Growth restriction parameters -------------------------------------------
mu = 0;
restrictVector = [];
restrictLengths = [];

restrictGrowth = false;

% Minimization parameters -------------------------------------------------
param = struct();
param.m = 20;
param.epsilon = 1e-5;
param.epsilon_rel = 1e-7;
param.past = 0;
param.delta = 0;
param.max_iterations = 10000;
param.max_linesearch = 20;
param.linesearch = 3;
param.linesearch_method = 1;
param.min_step = 1e-20;
param.max_step = 1e20;
param.ftol = 1e-4;
param.wolfe = 0.9;
param.iterDisp = 1;

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if(nargin<1), error('Please supply face connectivity list!'); end
if(nargin<2), error('Please supply vertex coordinates!'); end
if(nargin<3), error('Please supply target edge lengths!'); end

Nv = size(vertex, 1); % The number of vertices
% Nf = size(face, 1); % The number of faces

% Validate the input triangulation
validateattributes(vertex, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'} );
validateattributes(face, {'numeric'}, ...
    {'2d', 'ncols', 3, 'positive', 'finite', ...
    'integer', 'real', '<=', Nv} );

% Check the size of the edge list
tr = triangulation( face, vertex );
edgeList = tr.edges;
v1 = edgeList(:,1);
v2 = edgeList(:,2);

% Construct topological structure tools for the mesh
e1IDx = sort( [ face(:,3), face(:,2) ], 2 );
e2IDx = sort( [ face(:,1), face(:,3) ], 2 );
e3IDx = sort( [ face(:,2), face(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, edgeList, 'rows' );
[~, e2IDx] = ismember( e2IDx, edgeList, 'rows' );
[~, e3IDx] = ismember( e3IDx, edgeList, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

% Default curvature is plate-like (i.e., flat) ----------------------------
tarTheta = zeros( size( tarLength ) );

% Check the size of the target length list
sizee = size( v1 );
sizetl = size( tarLength );
if ( ~isequal( sizee, sizetl ) )
    tarLength = tarLength';
    sizetl = size(tarLength);
    if ( ~isequal( sizee, sizetl ) )
        error( 'The target length list is improperly sized');
    end
end

for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
    
    % Target geometry parameters ------------------------------------------
    if strcmpi(varargin{i}, 'TargetAngles')
        tarTheta = varargin{i+1};
    end
    
    % Target vertex correspondence parameters -----------------------------
    if strcmpi(varargin{i}, 'FixBoundary')
        fixBoundary = true;
        anyFixed = true;
        alpha = 1;
    end    
    if strcmpi(varargin{i}, 'TargetVertices')
        target_ID = varargin{i+1};
        anyFixed = true;
        alpha = 1;
    end
    if strcmpi(varargin{i}, 'TargetLocations')
        targetLocations = varargin{i+1};
        anyFixed = true;
        alpha = 1;
    end
    if strcmpi(varargin{i}, 'Alpha')
        alpha = varargin{i+1};
    end
    
    % Target volume correspondence parameters -----------------------------
    if strcmpi(varargin{i}, 'FixVolume')
        fixVolume = true;
        beta = 1;
    end
    if strcmpi(varargin{i}, 'TargetVolume')
        targetVolume = varargin{i+1};
        fixVolume = true;
        beta = 1;
    end
    if strcmpi(varargin{i}, 'Beta')
        beta = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'PhantomFaces')
        phantomFaces = varargin{i+1};
        usePhantom = true;
    end
    
    % Growth restriction parameters ---------------------------------------
    if strcmpi(varargin{i}, 'GrowthRestrictionField')
        restrictVector = varargin{i+1};
        restrictGrowth = true;
    end
    if strcmpi(varargin{i}, 'RestrictedLengths')
        restrictLengths = varargin{i+1};
        restrictGrowth = true;
    end
    if strcmpi(varargin{i}, 'Mu')
        mu = varargin{i+1};
    end
    
    % Material parameters -------------------------------------------------
    if strcmpi(varargin{i}, 'Thickness')
        h = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'Poisson')
        nu = varargin{i+1};
    end
    
    % Minimization parameters ---------------------------------------------
    if strcmpi(varargin{i}, 'MHistorySize')
        param.m = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'Epsilon')
        param.epsilon = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'EpsilonRel')
        param.epsilon_rel = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'Past')
        param.past = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'Delta')
        param.delta = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'MaxIterations')
        param.max_iterations = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'MaxLinesearch')
        param.max_linesearch = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'MinStep')
        param.min_step = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'MaxStep')
        param.max_step = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'FTol')
        param.ftol = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'WolfeCoeff')
        param.wolfe = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'IterDisplay')
        if( varargin{i+1} )
            param.iterDisp = 1;
        else
            param.iterDisp = 0;
        end
    end
    if strcmpi(varargin{i}, 'Linesearch')
        if strcmpi(varargin{i+1}, 'Armijo')
            param.linesearch = 1;
        elseif strcmpi(varargin{i+1}, 'Wolfe')
            param.linesearch = 2;
        elseif strcmpi(varargin{i+1}, 'StrongWolfe')
            param.linesearch = 3;
        else
            error( ['Invalid line search termination ' ...
                'condition supplied!'] );
        end
    end
    if strcmpi(varargin{i}, 'LinesearchMethod')
        if strcmpi(varargin{i+1}, 'Backtracking')
            param.linesearch_method = 1;
        elseif strcmpi(varargin{i+1}, 'Bracketing')
            param.linesearch_method = 2;
        elseif strcmpi(varargin{i+1}, 'NocedalWright')
            param.linesearch_method = 3;
        else
            error('Invalid line search method supplied!');
        end
    end
    
end

if ( (param.linesearch_method == 3) && (param.linesearch ~= 3) )
    error(['The line search method of Nocedal and Wright can ' ...
        'only be used with the Strong Wolfe termination conditions']);
end

% Process target geometry input -------------------------------------------

% Check that the target edge list satisfies the triangle inequality
tarL_F = tarLength(feIDx);
assert( all( (sum(tarL_F, 2) - 2 .* tarL_F) > 0 ), ...
    'Target edge length list does not satisfy the triangle inequality' );

% Check the size of the target angle list
sizeta = size( tarTheta );
if ( ~isequal( sizee, sizeta ) )
    tarTheta = tarTheta';
    sizeta = size(tarTheta);
    if ( ~isequal( sizee, sizeta ) )
        error( 'The target angle list is improperly sized');
    end
end

% Process growth restriction input ----------------------------------------

if restrictGrowth
    
    % Check that the value of mu is valid
    assert( mu > 0, ...
        'Mu must be positive if a growth restriction field is supplied' );
    
    % Check that all necessary inputs have been supplied
    assert( ~isempty(restrictVector) && ~isempty(restrictLengths), ...
        ['The restriction field and the restricted ' ...
        'lengths must both be supplied'] );
    
    % Check that the restriction field is valid
    validateattributes(restrictVector, {'numeric'}, ...
        {'2d', 'nrows', size(face,1), 'ncols', 3, 'finite', 'real'});
    
    % Check that the restriction lengths are valid
    validateattributes(restrictLengths, {'numeric'}, ...
        {'2d', 'nrows', size(face,1), 'ncols', 3, ...
        'positive', 'finite', 'real'});
    
    % Ensure that the restriction field is tangent to the surface
    fN = faceNormal(tr);
    restrictVector = restrictVector - ...
        repmat(dot(restrictVector, fN, 2), 1, 3) .* fN;
    
    % Normalize the restriction vector
    restrictVector = restrictVector ./ ...
        repmat( sqrt(sum(restrictVector.^2, 2)), 1, 3 );
    
    % Check that no input edge in the initial geometry already exceeds its
    % maximum length along the restriction field
    edgeVecX = vertex(edgeList(:,2), 1) - vertex(edgeList(:,1), 1);
    edgeVecY = vertex(edgeList(:,2), 2) - vertex(edgeList(:,1), 2);
    edgeVecZ = vertex(edgeList(:,2), 3) - vertex(edgeList(:,1), 3);
    
    edgeVecF = cat(3, edgeVecX(feIDx), edgeVecY(feIDx), ...
        edgeVecZ(feIDx) );
    
    projL = repmat(permute(restrictVector, [1 3 2]), 1, 3, 1);
    projL = abs(dot(projL, edgeVecF, 3));
    
    assert( all( projL(:) < restrictLengths(:) ), ...
        [ 'Some edges projected along the restriction field already ' ...
        'exceed their maximum allowed lengths' ] );
    
else
    
    % Prepare non-empty variables for C++ functionality check
    restrictVector = -ones(1,3);
    restrictLengths = -1;
    
end


% Process target volume correspondence input ------------------------------

if fixVolume
    
    % Check that the value of beta is valid
    assert( beta > 0, ...
        'Beta must be positive if a target volume is supplied' );
    
    if usePhantom
        
        validateattributes(phantomFaces, {'numeric'}, ...
            {'2d', 'ncols', 3, 'finite', 'positive', 'real', 'integer'});
        
        assert( isequal(unique(phantomFaces), unique(face)), ...
            'Invalid phantom face connectivity list' );
        
        % Check that the supplied surface has no boundaries
        assert( ...
            isempty(freeBoundary(triangulation(phantomFaces, vertex))), ...
            'Phantom triangulation must represent a closed surface' );
        
    else
        
        % Check that the supplied surface has no boundaries
        assert( isempty(freeBoundary(tr)), ...
            [ 'Only surfaces without boundary can be assigned ' ...
            'a target volume. Consider using a closed phantom ' ...
            'triangulation' ] );
        
    end
    
    % Calculate target volume if necessary
    if isempty(targetVolume)
        
        if usePhantom
            
            % The centroids of each face
            COM = cat( 3, vertex(phantomFaces(:,1), :), ...
                vertex(phantomFaces(:,2), :), ...
                vertex(phantomFaces(:,3), :) );
            COM = mean(COM, 3);
            
            % The area weighted face normal vectors
            ej = vertex(phantomFaces(:,1), :) - ...
                vertex(phantomFaces(:,3), :);
            ek = vertex(phantomFaces(:,2), :) - ...
                vertex(phantomFaces(:,1), :);
            n = cross(ej, ek, 2);
            
        else
            
            % The centroids of each face
            COM = cat( 3, vertex(face(:,1), :), ...
                vertex(face(:,2), :), vertex(face(:,3), :) );
            COM = mean(COM, 3);
            
            % The area weighted face normal vectors
            ej = vertex(face(:,1), :) - vertex(face(:,3), :);
            ek = vertex(face(:,2), :) - vertex(face(:,1), :);
            n = cross(ej, ek, 2);
            
        end
        
        targetVolume = sum( dot(COM, n, 2) ) ./ 6;
        
        if any(targetVolume < 0)
            
            targetVolume = -targetVolume;
            face = face(:, [2 1 3]);
            
            if usePhantom
                phantomFaces = phantomFaces(:, [2 1 3]);
            end
            
        elseif any(targetVolume == 0)
            
            error('Invalid user supplied mesh');
            
        end
        
    else
        
        % Validate user supplied target volume
        assert( targetVolume > 0, 'Target volume must be positive' );
        
    end
    
else
    
    % Prepare non-empty variables for C++ functionality check
    targetVolume = -1;
    phantomFaces = -ones(1, 3);
    
end


% Process target vertex correspondence input ------------------------------

% Check that the value of alpha is valid
if ( anyFixed && ( alpha <= 0 ) )
    error( 'Alpha must be positive if target vertices are supplied' );
end

% Check that the target ID list and the target location list are consistent
if ( size( target_ID, 1 ) ~= size( targetLocations, 1 ) )
    error( [ 'Size of target vertex ID list does not', ...
        'match the size of the target location list' ] );
end

% Check the size of the target location list
if ( ~isempty(targetLocations) && ( size( targetLocations, 2 ) ~= 3 ) )
    error( 'Target location list is improperly sized' );
end

% Construct fixed boundary condition if necessary
if ( fixBoundary )
    
    bdyID = unique( tr.freeBoundary, 'stable' );
    
    target_ID = [ target_ID; bdyID ];
    targetLocations = [ targetLocations; vertex( bdyID, : ) ];
    
end

% Check if target vertex indices exist
if ( anyFixed )
    
    if ( max(target_ID(:)) > size(vertex,1) )
        error('The target ID list contains an undefined vertex');
    elseif ( min(target_ID(:)) < 1 )
        error('The target ID list contains a vertex index smaller than 1');
    end
    
else
    
    % Prepare non-empty array for C++ functionality check
    target_ID = -1;
    targetLocations = zeros(1,3);
    
end


% Update vIDs to match the 0-indexing in C++
face = face-1; v1 = v1-1; v2 = v2-1;

if (usePhantom)
    phantomFaces = phantomFaces-1;
end

if (anyFixed) 
    target_ID = target_ID-1;
end

% -------------------------------------------------------------------------
% Run Elastic Minimization!
%--------------------------------------------------------------------------

newVertex = minimize_elastic_energy( int32(face), double(vertex), ...
    int32(v1), int32(v2), double(tarLength), double(tarTheta), ...
    param, double(h), double(nu), ...
    double(alpha), int32(target_ID), double(targetLocations), ...
    double(beta), double(targetVolume), ...
    usePhantom, int32(phantomFaces), ...
    double(mu), double(restrictVector), double(restrictLengths) );

newVertex = reshape(newVertex, [Nv 3]);

end

