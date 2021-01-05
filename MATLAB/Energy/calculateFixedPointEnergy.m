function Efp = calculateFixedPointEnergy(F, V, fixedIDx, fixedX, alpha)
%CALCULATEFIXEDPOINTENERGY Calculates the fixed point energy
%used in the Non-Euclidean Shell Simulator (NES) for a given configuration
%and set of fixed vertices and fixed vertex locations
%
%   INPUT PARAMETERS:
%
%       F:          #Fx3 face connectivity list
%       V:          #Vx3 3D vertex coordinate list
%       fixedIDx:   #Px1 list of fixed vertex IDs
%       fixedX:     #Px3 list of fixed vertex target coordinates
%       alpha:      The scalar coefficient of the fixed point energy
%
%   OUTPUT PARAMETERS:
%
%       Efp:     The fixed point energy
%
%   by Dillon Cislo 2021/01/02

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinate list'); end
if (nargin < 3), error('Please supply fixed vertex ID list'); end
if (nargin < 4), error('Please supply target coordinates'); end
if (nargin < 5), alpha = 1; end

validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan'} );
validateattributes( fixedIDx, {'numeric'}, ...
    {'vector', 'integer', 'positive'} );
validateattributes( fixedX, {'numeric'}, ...
    {'2d', 'ncols', 3, 'nrows', numel(fixedIDx), 'finite', 'nonnan'} );
validateattributes( alpha, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan', '>=', 0});

if (size(fixedIDx,2) ~= 1), fixedIDx = fixedIDx.'; end

%--------------------------------------------------------------------------
% Calculate the Fixed Point Energy
%--------------------------------------------------------------------------

Efp = alpha * sum(sum((V(fixedIDx, :) - fixedX).^2, 2));

end

