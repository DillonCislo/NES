function [Vol, outwardNormal] = calculateMeshVolume(F, V)
%CALCULATEMESHVOLUME Calculates the volume enclosed by a boundaryless mesh
%triangulation
%
%   INPUT PARAMETERS:
%       - F:                #Fx3 face connectivity list
%       - V:                #VxD vertex coordinate list
%
%   OUTPUT PARAMETERS:
%       - Vol:              The volume enclosed by the mesh
%       - outwardNormal:   	boolean indicating whether the face normals
%                           calculated from the input face connectivity
%                           list are outward pointing
%
%   by Dillon Cislo 03/05/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});

TR = triangulation(F, V);

assert( isempty(TR.freeBoundary), ...
    'Enclosed volumes are not well defined for surfaces with boundaries' );

%--------------------------------------------------------------------------
% CALCULATE VOLUME
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
if any(Vol < 0)
    Vol = -Vol;
    outwardNormal = false;
else
    outwardNormal = true;
end
    

end

