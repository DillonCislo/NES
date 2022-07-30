function efIDx = edgeFaceCorrespondence(F)
%EDGEFACECORRESPONDENCE Determines the faces attached to each edge in a
%mesh triangulation. Ordering or the faces in the output reflects the
%ordering of faces used to calculate the bending angles in NES
%
%   INPUT PARAMETERS:
%
%       - F:         #Fx3 face connectivity list
%
%   OUTPUT PARAMETERS:
%
%       - efIDx:        #Ex2 edge-face correspondence list. efIDx(i,j)
%                       refers to the jth face attached to the ith edge
%
%   by Dillon Cislo 02/01/2022

validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', 'finite'});

% Generate a set of phantom vertices
V = zeros(max(F(:)), 2);

TR = triangulation(F,V);
E = TR.edges;
v1 = E(:,1);
v2 = E(:,2);

efIDx = edge_face_correspondence(F, V, v1, v2);

end