
## adaptiveRemesh.m

Performs the adaptive remesh. qOcta can either be ```[]``` for no initialization, or a coarse qOcta that corresponds to the model. saveIterates is a boolean for whether or not in-depth info is needed, and addInfo is present in case one wishes to append to an existing info struct. It returns the new info struct, the new qOcta, and the new model.

## edgeSplitXT.m

Performs the edge split, given a ```model```, its associated ```costs```, and a ```qOcta```. ```model``` is expected to be a valid model, as described in README.md. ```costs``` is expected to have length corresponding to the number of fine control points (both vertices and edges), although the points at the vertex will always have a cost of zero, such that the ```i```th entry corresponds to the control point with index ```i``` in the model. ```qOcta``` should be the coarse value of q.

## euclideanCost.m

Gets the cost by taking Euclidean distance.