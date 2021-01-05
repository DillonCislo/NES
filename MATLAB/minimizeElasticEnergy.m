function [newVertex] = minimizeElasticEnergy(face, vertex, tarLength, ...
    varargin )
%MINIMIZEELASTICENERGY this function minimizes the elastic energy of a
%non-Euclidean shell with a given initial configuration and a specified
%intrinsic geometry
%   Detailed explanation goes here

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
param.max_iterations = 20000;
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

% Check the size of the face connectivity list
sizef = size(face);
if ((sizef(2)~=3)||(length(sizef)~=2))
    error('The face list is improperly sized');
end

% Check the size of the vertex list
sizev = size(vertex);
if ((sizev(2)~=3)||(length(sizev)~=2))
    error('The vertex list is improperly sized');
end

% Check if vertex indices exist
if ( max(face(:)) > sizev(1) )
    error('The face list contains an undefined vertex');
elseif ( min(face(:)) < 1 )
    error('The face list contains a vertex index smaller than 1');
end

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
    
    % if any( projL(:) > restrictLengths(:) )
    %     
    %     warning( [ 'Some edges projected along the restriction ' ...
    %         'field already exceed their maximum allowed lengths' ] );
    %     
    % end
    
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
        
        if (targetVolume < 0)
            
            targetVolume = -targetVolume;
            face = face(:, [2 1 3]);
            
            if usePhantom
                phantomFaces = phantomFaces(:, [2 1 3]);
            end
            
        elseif (targetVolume == 0)
            
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

newVertex = reshape(newVertex, sizev);


end

