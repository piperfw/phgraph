# phgraph
Graph algorithms implemented using the Cytoscape.js library [https://js.cytoscape.org](https://js.cytoscape.org)

## Using the library
- Include `phgraph.js` in a `<script>` tag after loading the Cytoscape.js library
- Instantiate a cytoscape object using `phCytoscape()`:
```js
let cy = phCytoscape({ /* options */ });
``` 
where `{ /* options */ }` are the same set of options you would pass to `cytoscape()` to initialise a core cytoscape
object. In particular the `container` option must be set in order for a graph to be displayed.
For use of the core Cytoscape.js library, refer to [https://js.cytoscape.org](https://js.cytoscape.org).
- Functions defined in `phgraph.js` may be accessed via the resulting object, such as in
```js
let collection = cy.phMakeGraph(degreeSequence=[1,1,2]);
```
which creates a graph with degree sequence `[1,1,2]`. These functions may also be accessed directly, typically taking
a cytoscape instance or collection as their first argument.

## Available functions
`phgraph.js` implements a set of graph theory algorithms in addition to those provided by Cytoscape.js, as well as
functions to test and build graphs of certain types. These include:

#### Testing
- `cy.phHasGraphRealisation(degreeSequence)` - Test whether a degree sequence (array) is graphical
- `cy.phHasBipartiteRealisation(degreeSequence)` - Test whether a degree sequence (array) has a bipartite realisation
- `cy.phHasTreeRealisation(degreeSequence)` - Test whether a degree sequence (array) has a tree realisation
- `cy.phTreesIsomorphic(collection1, collection2)` - Test two collections describe isomorphic trees

#### Construction
- `cy.phMakeGraph(degreeSequence)` - Add a graph to `cy` with a given degree sequence (array)
- `cy.phMakeBipartite(degreeSeq1, degreeSeq2)` - Add a bipartite graph to `cy` with a given degree sequence (two arrays)
- `cy.phMakeTree(degreeSequence)` - Add a tree to `cy` with a given degree sequence (array)
- `cy.phMakeComplement()` - Replace the current graph in `cy` with its complement
- `cy.phDegreeSequence()` - Return the degree sequence (array) for the graph currently in `cy`

#### Other
- `cy.phListFunctions()` - Return array of functions bound to `cy`
- `cy.phEnumerateTrees(degreeSequence)` - Attempt to enumerate all trees with a given degree sequence, up to isomorphism


Comments in `phgraph.js` provide brief documentation and a description of optional parameters for these functions.
For example usage see `index.html`.