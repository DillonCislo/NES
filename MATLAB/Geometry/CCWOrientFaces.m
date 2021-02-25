function [Fnew, isCCW, isCW] = CCWOrientFaces(F, V)
%CCWORIENTFACES Re-orders the face connectivity list of an input mesh
%triangulation so that all faces are counter-clockwise ordered in the plane
%of the face
%
%   INPUT PARAMETERS:
%
%       - F:        #Fx3 face connectivity list
%       - V:        #VxD vertex coordinate list
%
%   OUTPUT PARAMETERS:
%
%       - Fnew:     #Fx3 re-ordred face connectivity list
%       - isCCW:    Is true if all faces of the input connectivity list
%                   are already CCW ordered
%       - isCW:     Is true if all faces of the input connectivity list
%                   are already CW ordered
%
% by Dillon Cislo 02/24/2021 

% Validate Inputs ---------------------------------------------------------

validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});

% Determine Current Face Orientation --------------------------------------

% The extrinsic edge vector for the edge opposite vertex k
Ek = V(F(:,2), :) - V(F(:,1), :);
Lk = sqrt( sum( Ek.^2, 2 ) ); % The length of the edge

% The extrinsic edge vector for the edge opposite vertex j
Ej = V(F(:,3), :) - V(F(:,1), :);
Lj = sqrt( sum( Ej.^2, 2 ) ); % The length of the edge

% The length of the edge vector opposite vertex i
Li = V(F(:,3), :) - V(F(:,2), :);
Li = sqrt( sum( Li.^2, 2 ) );

% Functions of the internal angle associate to vertex i
cosAng = ( Lj.^2 + Lk.^2 - Li.^2 ) ./ ( 2 .* Lj .* Lk );
sinAng = sin(acos(cosAng));

% The intrinsic edge vector for the edge opposite vertex k
ek = [ Lk, zeros(size(F,1), 2) ];

% The intrinsic edge vector for the edge opposite vertex j
ej = Lj .* [ cosAng, sinAng, zeros(size(F,1), 1) ];

% Faces are CCW ordered if the sign of the 3rd component of the
% cross product of the intrinsic edges are positive
faceOrder = cross(ek, ej, 2);
faceOrder = sign( faceOrder(:, 3) );

% Update Face Orientation -------------------------------------------------

isCCW = all(faceOrder > 0);
isCW = all(faceOrder < 0);

Fnew = F;
for f = 1:size(F,1)
    
    if faceOrder(f) < 0
        Fnew(f,:) = fliplr(Fnew(f,:));
    end
    
end

end

