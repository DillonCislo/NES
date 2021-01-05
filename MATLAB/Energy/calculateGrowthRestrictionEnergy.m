function [Erg, projL, isValid] = calculateGrowthRestrictionEnergy(F, V, ...
    growthVec, maxProjL, mu)
%CALCULATEFIXEDPOINTENERGY Calculates the grwoth restriction energy
%used in the Non-Euclidean Shell Simulator (NES) for a given configuration
%and a set of maximum projected edge lengths
%
%   INPUT PARAMETERS:
%
%       F:              #Fx3 face connectivity list
%       V:              #Vx3 current 3D vertex coordinate list
%       growthVec:      #Fx3 tangent unit vector field mesh faces
%       maxProjL:       #Fx3 maximum projected length of each edge along
%                       the growth unit vector field in each face
%       mu:             The inverse scalar coefficient of the growth
%                       restriction energy
%
%   OUTPUT PARAMETERS:
%
%       Erg:            The growth restriction energy
%       projL:          The projected edge lengths
%       isValid:        If true, all projected edge lengths are less than
%                       their maximum allowed length
%
%   by Dillon Cislo 2021/01/02

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinate list'); end
if (nargin < 3), error('Please supply growth restriction field'); end
if (nargin < 4), error('Please supply max projected edge lengths'); end
if (nargin < 5), mu = 1; end

validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan'} );
validateattributes( growthVec, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'nrows', size(F,1)} );
validateattributes( maxProjL, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'positive', ...
    'nonnan', 'nrows', size(F,1)} );
validateattributes( mu, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan', '>=', 0});

TR = triangulation(F, V);
E = TR.edges;

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, E, 'rows' );
[~, e2IDx] = ismember( e2IDx, E, 'rows' );
[~, e3IDx] = ismember( e3IDx, E, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

%--------------------------------------------------------------------------
% Calculate Growth Restriction Energy
%--------------------------------------------------------------------------

% Cartesian components of directed edge vectors
edgeVecX = V(E(:,2), 1) - V(E(:,1), 1);
edgeVecY = V(E(:,2), 2) - V(E(:,1), 2);
edgeVecZ = V(E(:,2), 3) - V(E(:,1), 3);

% Edge vectors defined on faces
% Each (#Fx3) page corresponds to a single Cartesian component
edgeVecF = cat(3, edgeVecX(feIDx), edgeVecY(feIDx), edgeVecZ(feIDx) );

% Calculate the projected edge lengths
projL = repmat(permute(growthVec, [1 3 2]), 1, 3, 1);
projL = abs(dot(projL, edgeVecF, 3));

isValid = all(projL(:) < maxProjL(:));

Erg = -sum(sum(log(maxProjL - projL), 2), 1) ./ mu;

end

