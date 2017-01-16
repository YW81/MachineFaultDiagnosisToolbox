function [inv varargout] = myPInv(mat)
%--------------------------------------------------------------------------
% Syntax:       [inv isSingular] = myPInv(mat);
%               inv = myPInv(mat);
%
% Inputs:       mat is an arbitrary M x N matrix
%
% Outputs:      inv is the N x M Moore-Penrose pseudoinverse of input mat.
%               When mat is square and invertible, inv = mat^(-1)
%
%               isSingular = 'true' for all nonsquare input mat and for
%               singular square input mat
%               isSingular = 'false' for all square invertible input mat
%              
% Description:  Computes the N x M Moore-Penrose pseudoinverse of the
%               M x N input mat using the numerically stable compact
%               singular value decomposition.
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         June 28, 2012
%--------------------------------------------------------------------------

[U S V] = mySVD(mat,'compact');
r = size(S,1);
if (r < max(size(mat)))
    isSingular = 'true';
else
    isSingular = 'false';
end

if r > 0
    inv = V * diag(ones(r,1) ./ diag(S)) * U';
else
    inv = zeros(size(mat,2),size(mat,1));
end

switch nargout
    case 1
        varargout = {};
    case 2
        varargout{1} = isSingular;
end
