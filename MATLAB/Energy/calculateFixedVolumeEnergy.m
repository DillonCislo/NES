function Efv = calculateFixedVolumeEnergy(F, V, tarVol, beta)
%CALCULATEFIXEDPOINTENERGY Calculates the fixed volume energy
%used in the Non-Euclidean Shell Simulator (NES) for a given configuration
%and target volume
%
%   INPUT PARAMETERS:
%
%       F:          #Fx3 face connectivity list
%       V:          #Vx3 3D vertex coordinate list
%       tarVol:     The target volume of the enclosed surface
%       beta:       The scalar coefficient of the fixed volume energy
%
%   OUTPUT PARAMETERS:
%
%       Efv:     The fixed volume energy
%
%   by Dillon Cislo 2021/01/02

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinate list'); end
if (nargin < 3), error('Please target volume'); end
if (nargin < 4), beta = 1; end

validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan'} );
validateattributes( tarVol, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan', 'positive'});
validateattributes( beta, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan', '>=', 0});

TR = triangulation(F, V);

if ~isempty(TR.freeBoundary)
    
    try
        
        FF = fill_holes(V, F);
        warning('Filling holes in mesh to calculate volume');
        F = FF;
        
    catch
        
        error(['Enclosed volumes are not well defined ' ...
            'for surfaces with boundaries']);
        
    end
    
end

%--------------------------------------------------------------------------
% Calculate Fixed Volume Energy
%--------------------------------------------------------------------------

% The centroids of each face
COM = mean( cat( 3, V(F(:,1), :), V(F(:,2), :), V(F(:,3), :) ), 3 );

% The area weighted face normal vectors
ej = V(F(:,1), :) - V(F(:,3), :);
ek = V(F(:,2), :) - V(F(:,1), :);
n = cross(ej, ek, 2);

Vol = dot(COM, n, 2);

if (numel(unique(sign(Vol))) ~= 1)
    warning('Mesh faces do not appear to be consistently ordered');
end

% The enclosed volume
Vol = sum(Vol) ./ 6;

% Correct for inward facing normals
if any(Vol < 0), Vol = -Vol; end

Efv = beta .* (Vol - tarVol).^2;

end

