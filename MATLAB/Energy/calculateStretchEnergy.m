function [Es, EsF] = calculateStretchEnergy(F, V, tarL, nu)
%CALCULATESTRETCHENERGY Calculates the (re-scaled) stretching energy used
%in the Non-Euclidean Shell Simulator (NES) for a given configuration and
%set of target edge lengths
%
%   INPUT PARAMETERS:
%
%       F:      #Fx3 face connectivity list
%       V:      #Vx3 3D vertex coordinate list
%       tarL:   #Ex1 list of target edge lengths
%       nu:     The 2D Poisson ratio
%
%   OUTPUT PARAMETERS:
%
%       Es:     The stretch energy
%       EsF:    #Fx1 stretch energy on each face
%
%   by Dillon Cislo 2021/01/02

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinate list'); end
if (nargin < 3), error('Please supply target edge list'); end
if (nargin < 4), nu = 1/3; end

validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan'} );
validateattributes( tarL, {'numeric'}, ...
    {'vector', 'finite', 'nonnan'} );
validateattributes( nu, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan'});

TR = triangulation(F,V);
E = sort(TR.edges, 2);

if (size(tarL,2) ~= 1), tarL = tarL.'; end
if (numel(tarL) ~= numel(E(:,1)))
    error('Target edge length list is improperly sized');
end

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
% Calculate the Stretch Energy
%--------------------------------------------------------------------------

% Target edge lengths defined on mesh faces
tarL_F = tarL(feIDx);

% Calculate target face areas
tarA = sum(tarL_F, 2) ./ 2; % the triangle semi-perimeters
tarA = sqrt( tarA .* (tarA - tarL_F(:,1)) .* ...
    (tarA - tarL_F(:,2)) .* (tarA - tarL_F(:,3)) );

% Physical mesh edge lengths
L = V(E(:,2), :) - V(E(:,1), :);
L = sqrt(sum(L.^2, 2));
L_F = L(feIDx); % defined on mesh faces

% Calculate the cyclic strains on edges in mesh faces
cycStr_F = L_F.^2 - tarL_F.^2;
cycStr_F = cycStr_F - ...
    circshift(cycStr_F, 1, 2) - circshift(cycStr_F, 2, 2);

% Calculate the trace of the strain tensor on each face -------------------
TrStr = cycStr_F .* tarL_F.^2 ./ repmat(tarA.^2, 1, 3);
TrStr = -sum(TrStr, 2) ./ 8;

% Calculate the trace of the squared strain tensor on each face -----------
prodStr = [ cycStr_F .* repmat(cycStr_F(:,1), 1, 3), ...
    cycStr_F .* repmat(cycStr_F(:,2), 1, 3), ...
    cycStr_F .* repmat(cycStr_F(:,3), 1, 3) ];

dotProd = tarL_F.^2;
dotProd = dotProd - circshift(dotProd, 1, 2) - circshift(dotProd, 2, 2);
dotProd = dotProd ./ 2;
dotProd = [ tarL_F(:,1).^2, dotProd(:,3), dotProd(:,2), ...
    dotProd(:,3), tarL_F(:,2).^2, dotProd(:,1), ...
    dotProd(:,2), dotProd(:,1), tarL_F(:,3).^2 ];

TrStr2 = sum( prodStr .* dotProd.^2, 2 ) ./ (64 .* tarA.^4);

% Combine to form the stretch energy
EsF = tarA .* ( nu .* TrStr.^2 + (1-nu) .* TrStr2 ) ./ 8;
Es = sum(EsF);

% Es = sum(tarA .* ( nu .* TrStr.^2 + (1-nu) .* TrStr2 ));
% Es = Es ./ 8;

end

