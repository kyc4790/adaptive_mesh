# deriv/

## Format

The order of derivatives with respect to ```q``` is ```q[1,1], q[1,2], ..., q[1,9], q[2,1], ...```.

The order of the second derivatives will be in the form of ```(q[1,1], q[1,1]), (q[1,1],q[1,2]), (q[1,1],q[1,3])..., (q[1,2], q[1,1], ...)```.

The order of constraints is as follows:
- The first ```numBdry``` constraints correspond to ```|x|^2 - 1 = 0```, in the order of ```mesh.bdryIdx```
- The next ```7 * numBdry``` constraints correspond to ```W_n * x  = v```, such that constraints 1 through 7 are the constraints on the first boundary point, constraints 8-14 are constraints on the second, etc
- The last ```15 * numInt``` constraints correspond to ```x^T P x = 0```, such that all constraints associated with ```P_1``` are in the first ```numInt``` constraints, and so on, such that each set of ```numInt``` constraints follows the order of ```mesh.intIdx```

## Constraints

Both constraint functions take in 

### getSensitivityConstriaint.m

Uses the method presented in the new Overleaf document to compute costs. Currently does not work since the matrix C_y is not invertible.

### getDerivConstraint.m

Currently uses the envelope theorem to calculate the costs. It also contains commented out code corresponding to the second derivative constraint in the original Overleaf document.
It is called in AdaptiveRemesh.m to compute the costs.

## Tests

### finite_differences_test/

This contains all of the code for the finite differences tests.

### checkEnvelopeTheorem.m

This file prints out the dot product of the costs with the results from the envelope theorem.

## Helper functions

### buildIHh.m
Extends a solution from coarse to fine.

### getAngularMomentumMatrices.m
Returns the angular momentum matrices.

### OctaMat.mat
Contains all of the matrices for the octahedral constraint

### rowspaces.m
Computes the row space, null space, and dimensions of a given array.

