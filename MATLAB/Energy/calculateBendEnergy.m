function [Eb, EbF] = calculateBendEnergy(F, V, tarL, tarAng, nu, h)
%CALCULATEBENDENERGY Calculates the (re-scaled) bending energy used
%in the Non-Euclidean Shell Simulator (NES) for a given configuration and
%set of target edge lengths and hinge angles
%
%   INPUT PARAMETERS:
%
%       F:          #Fx3 face connectivity list
%       V:          #Vx3 3D vertex coordinate list
%       tarL:       #Ex1 list of target edge lengths
%       tarAng:     #Ex1 list of target edge hinge angles
%       nu:         The 2D Poisson ratio
%       h:          The thickness of the elastic body
%
%   OUTPUT PARAMETERS:
%
%       Eb:     The bend energy
%       EbF:    #Fx1 bend energy on each face
%
%   by Dillon Cislo 2021/01/02

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinate list'); end
if (nargin < 3), error('Please supply target edge length list'); end
if (nargin < 4), error('Please supply target hinge angle list'); end
if (nargin < 5), nu = 1/3; end
if (nargin < 6), h = 0.01; end

validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan'} );
validateattributes( tarL, {'numeric'}, ...
    {'vector', 'finite', 'nonnan'} );
validateattributes( tarAng, {'numeric'}, ...
    {'vector', 'finite', 'nonnan'} );
validateattributes( nu, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan'});
validateattributes( h, {'numeric'}, ...
    {'scalar', 'finite','nonnan', 'positive'});

TR = triangulation(F,V);
E = TR.edges;
v1 = E(:,1);
v2 = E(:,2);

if (size(tarL,2) ~= 1), tarL = tarL.'; end
if (numel(tarL) ~= numel(E(:,1)))
    error('Target edge length list is improperly sized');
end

if (size(tarAng,2) ~= 1), tarAng = tarAng.'; end
if (numel(tarAng) ~= numel(E(:,1)))
    error('Target edge length list is improperly sized');
end

% Construct edge-face correspondence tool ---------------------------------
efIDx = TR.edgeAttachments(E);
efIDx = cellfun(@(x) repmat(x, 1, 1+mod(numel(x),2)), efIDx, 'Uni', false);
efIDx = cell2mat(efIDx);

bdyEdge = ((efIDx(:,1) - efIDx(:,2)) == 0);

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
% Calculate Bend Energy
%--------------------------------------------------------------------------

% Target edge lengths defined on mesh faces
tarL_F = tarL(feIDx);

% Target hinge angle defined on mesh faces
tarAng_F = tarAng(feIDx);

% The hinge function defined on mesh faces
tarPhi_F = 2 .* tan( tarAng_F ./ 2 );

% Calculate target face areas
tarA = sum(tarL_F, 2) ./ 2; % the triangle semi-perimeters
tarA = sqrt( tarA .* (tarA - tarL_F(:,1)) .* ...
    (tarA - tarL_F(:,2)) .* (tarA - tarL_F(:,3)) );

% Calculate physical mesh hinge functions----------------------------------
% phi(theta) = 2 * tan( theta / 2 )

% OLD WAY - PRODUCES INCORRECT ANGLE SIGN FOR SOME EDGES ******************
% fN = TR.faceNormal;
% N1 = fN(efIDx(:,1), :);
% N2 = fN(efIDx(:,2), :);
% 
% crossN = cross(N1, N2, 2);
% crossHat = crossN ./ sqrt(sum(crossN.^2, 2));
% 
% crossN(bdyEdge, :) = zeros(sum(bdyEdge), 3);
% crossHat(bdyEdge, :) = zeros(sum(bdyEdge), 3);
% 
% phi = 2 .* dot(crossN, crossHat, 2) ./ (1 + dot(N1, N2, 2));
% phi_F = phi(feIDx);
%**************************************************************************

% Edge bending angles
% theta = calculate_edge_angles( int32(F-1), double(V), ...
%     int32(v1-1), int32(v2-1) );
theta = calculate_edge_angles( F, V, v1, v2 );

phi = 2 * tan( theta / 2 );
phi_F = phi(feIDx);

% Calculate bend strains on mesh faces
Phi_F = phi_F - tarPhi_F;

% Calculate the trace of the bend moment tensor on each face --------------
TrBM = sum( Phi_F .* tarL_F, 2 ) ./ ( 2 .* tarA );

% Calcualte the trace of the squared bend moment tensor on each face ------
prodPhi = [ Phi_F .* repmat(Phi_F(:,1), 1, 3), ...
    Phi_F .* repmat(Phi_F(:,2), 1, 3), ...
    Phi_F .* repmat(Phi_F(:,3), 1, 3) ];

dotProd = tarL_F.^2;
dotProd = dotProd - circshift(dotProd, 1, 2) - circshift(dotProd, 2, 2);
dotProd = dotProd ./ 2;
dotProd = [ tarL_F(:,1).^2, dotProd(:,3), dotProd(:,2), ...
    dotProd(:,3), tarL_F(:,2).^2, dotProd(:,1), ...
    dotProd(:,2), dotProd(:,1), tarL_F(:,3).^2 ];

prodL = [ tarL_F .* repmat(tarL_F(:,1), 1, 3), ...
    tarL_F .* repmat(tarL_F(:,2), 1, 3), ...
    tarL_F .* repmat(tarL_F(:,3), 1, 3) ];

TrBM2 = sum( prodPhi .* dotProd.^2 ./ prodL, 2 ) ./ (4 .* tarA.^2);

% Combine to form the bend energy
EbF = h.^2 .* tarA .* ( nu .* TrBM.^2 + (1-nu) .* TrBM2 ) ./ 24;
Eb = sum(EbF);

% Eb = sum(tarA .* ( nu .* TrBM.^2 + (1-nu) .* TrBM2 ));
% Eb = h.^2 .* Eb ./ 24;

end

