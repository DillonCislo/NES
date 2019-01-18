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

% Minimization parameters -------------------------------------------------
param = struct();
param.m = 20;
param.epsilon = 1e-7;
param.past = 0;
param.delta = 0;
param.max_iterations = 20000;
param.max_linesearch = 20;
param.linesearch = 3;
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
if ( max(face(:)) > size(vertex,1) )
    error('The face list contains an undefined vertex');
elseif ( min(face(:)) < 1 )
    error('The face list contains a vertex index smaller than 1');
end

% Check the size of the edge list
tr = triangulation( face, vertex );
edgeList = tr.edges;
v1 = edgeList(:,1);
v2 = edgeList(:,2);

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
    if ~isempty(regexp(varargin{i},'^[Tt]arget[Aa]ngles', 'match'))
        tarTheta = varargin{i+1};
    end
    
    % Target vertex correspondence parameters -----------------------------
    if ~isempty(regexp(varargin{i}, '^[Ff]ix[Bb]oundary', 'match'))
        fixBoundary = true;
        anyFixed = true;
        alpha = 1;
    end    
    if ~isempty(regexp(varargin{i}, '^[Tt]arget[Vv]ertices', 'match'))
        target_ID = varargin{i+1};
        anyFixed = true;
        alpha = 1;
    end
    if ~isempty(regexp(varargin{i}, '^[Tt]arget[Ll]ocations', 'match'))
        targetLocations = varargin{i+1};
        anyFixed = true;
        alpha = 1;
    end
    if ~isempty(regexp(varargin{i}, '^[Aa]lpha', 'match'))
        alpha = varargin{i+1};
    end
    
    % Material parameters -------------------------------------------------
    if ~isempty(regexp(varargin{i}, '^[Tt]hickness', 'match'))
        h = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]oisson', 'match'))
        nu = varargin{i+1};
    end
    
    % Minimization parameters ---------------------------------------------
    if ~isempty(regexp(varargin{i}, '^[Mm][Hh]istory[Ss]ize','match'))
        param.m = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ee]psilon','match'))
        param.epsilon = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]ast','match'))
        param.past = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Dd]elta','match'))
        param.delta = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Ii]terations', 'match'))
        param.max_iterations = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Ll]inesearch', 'match'))
        param.max_linesearch = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]in[Ss]tep', 'match'))
        param.min_step = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Ss]tep', 'match'))
        param.max_step = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ff][Tt]ol', 'match'))
        param.ftol = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ww]olfe[Cc]oeff', 'match'))
        param.wolfe = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ii]ter[Dd]isplay', 'match'))
        if( varargin{i+1} )
            param.iterDisp = 1;
        else
            param.iterDisp = 0;
        end
    end
    if ~isempty(regexp(varargin{i}, '^[Ll]inesearch', 'match'))
        switch varargin{i+1}
            case 'Armijo'
                param.linesearch = 1;
            case 'Wolfe'
                param.linesearch = 2;
            case 'StrongWolfe'
                param.linesearch = 3;
            otherwise
                error('Invalid line search method supplied!' );
        end
    end
    
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

% Process target vertex correspondence input ------------------------------

% Check that the value of alpha is consistent
if ( anyFixed && ( alpha == 0 ) )
    error( 'Alpha cannot equal zero if target vertices are supplied!' );
end

% Check that the target ID list and the target location list are consistent
if ( size( target_ID, 1 ) ~= size( targetLocations, 1 ) )
    error( [ 'Size of target vertex ID list does not',
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

if (anyFixed) 
    target_ID = target_ID-1;
end

% -------------------------------------------------------------------------
% Run Elastic Minimization!
%--------------------------------------------------------------------------

newVertex = minimize_elastic_energy( int32(face), double(vertex), ...
    int32(v1), int32(v2), double(tarLength), double(tarTheta), ...
    param, double(h), double(nu), ...
    double(alpha), int32(target_ID), double(targetLocations) );


end

