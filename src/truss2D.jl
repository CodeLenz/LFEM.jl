#
# Local stiffness matrix
#
function K_truss2D(E,A,L)
       (E*A/L)*[ 1.0 0.0 -1.0 0.0 ;
                 0.0 0.0  0.0 0.0 ; 
                -1.0 0.0  1.0 0.0 ;
                 0.0 0.0  0.0 0.0 ]
end
        