function [L, T] = calculateEdgeLengthsAndAngles(F, V)
%CALCULATE_EDGE_LENGTHS_AND_ANGLES Calculates the edge lengths and bending
%angles of a mesh triangulation.  Here the bending angles are defined as
%the angle between the normal vectors of the faces sharing an edge.
%
%   INPUT PARAMETERS:
%       - F:        #Fx3 face connectivity list
%       - V:        #VxD vertex coordinate list
%
%   OUTPUT PARAMETERS:
%       - L:        #Ex1 edge length list
%       - T:        #Ex1 edge bending angle list
%
%   by Dillon Cislo 03/05/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});

TR = triangulation(F, V);
E = TR.edges;
v1 = E(:,1);
v2 = E(:,2);

%--------------------------------------------------------------------------
% CALCULATE EDGE LENGTHS
%--------------------------------------------------------------------------

L = V(E(:,2), :) - V(E(:,1), :);
L = sqrt(sum(L.^2, 2));

%--------------------------------------------------------------------------
% CALCULATE BENDING ANGLES
%--------------------------------------------------------------------------

% OLD METHOD - PRODUCES INCORRECT ANGLE SIGN FOR SOME EDGES ---------------
% % #Ex2 array of fIDs of the faces attached to a particular edge.
% % If an edge is a border edge (i.e., only attached to a single face),
% % then that fID is listed twice for dimensional consistency
% resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
% edgeFace = edgeAttachments( TR, E );
% edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
%     'UniformOutput', false ) );
% 
% bdyEdge = ((edgeFace(:,1) - edgeFace(:,2)) == 0);
% 
% % Face unit normal vectors
% N = TR.faceNormal;
% N1 = N(edgeFace(:,1), :);
% N2 = N(edgeFace(:,2), :);
% 
% crossN = cross(N1, N2, 2);
% crossHat = crossN ./ sqrt(sum(crossN.^2, 2));
% 
% crossN(bdyEdge, :) = zeros(sum(bdyEdge), 3);
% crossHat(bdyEdge, :) = zeros(sum(bdyEdge), 3);
% 
% OldT = 2 .* atan2( dot(crossN, crossHat, 2), 1 + dot(N1, N2, 2) );
%--------------------------------------------------------------------------

% T = calculate_edge_angles( int32(F-1), double(V), ...
%     int32(v1-1), int32(v2-1) );

T = calculate_edge_angles( F, double(V), v1, v2 );

end

