% DPvar Toolbox
% Version 1.00, 9 September 2021.
%
% Creating dpvar objects
%    DPVAR         - Construct a dpvar
%
% Dpvar functions:
%    COMBINE       - Combine terms of dpvars
%    COMMON_BASIS  - Represent two dpvars using a shared basis
%    COMPRESS      - Remove terms with zero coefficients
%    DIFF          - Element-by-element differentiation of a dpvar
%    INT           - Element-by-element integration of a dpvar
%    JACOBIAN      - Compute Jacobian matrix of a vector dpvar
%    SUBS          - Symbolic substitution of independent variables
%    VARSWAP       - Symbolic swap of independent variable names
%
% Dpvar characteristics:
%    ISEMPTY       - True for empty dpvars
%    SIZE          - Size of a dpvar matrix
%    LENGTH        - Length of a dpvar matrix
%    GET           - Extract properties of a dpvar
%    SUBSREF       - Extract elements of a dpvar
%    SUBSASGN      - Assign elements of a dpvar matrix
%
% Conversions:
%    BSHAPE        - Convert coefficient matrix in multipoly format to dpvar
%                    format
%    DPVAR2POLY    - Convert object from multipoly toolbox to dpvar
%    POLY2DPVAR    - Convert object from dpvar to multipoly
%    DOUBLE        - Convert constant dpvar to a double
%
% Overloaded arithmetic operations:
%    PLUS, +       - Add dpvars
%    MINUS, -      - Subtract dpvars
%    MTIMES, *     - Multiply dparss
%    RDIVIDE, ./   - Elementwise division of dpvar by double
%    HORZCAT, [,]  - Horizontal concatentation of dpvars
%    VERTCAT, [;]  - Vertical concatentation of dpvars
%    BLKDIAG       - Block diagonal concatenation of dpvar matrices
%    CTRANSPOSE, ' - Non-conjugate transpose of a dpvar
%    TRANSPOSE, .' - Non-conjugate transpose of a dpvar
%    UPLUS         - Unary plus of a dpvar
%    UMINUS        - Unary minus of a dpvar
%

