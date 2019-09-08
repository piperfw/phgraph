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
- `js phListFunctions` - Return array of functions bound to `cy`
- 


Comments in `phgraph.js` provide minimal documentation for these functions; for example usage see `index.html`.