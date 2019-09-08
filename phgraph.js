/*
 * Wrapper to create a cytoscape instance with a set of additional methods defined in this document
 * @param {object literal} options - Options passed to cytoscape(). See http://js.cytoscape.org/#getting-started/initialisation 
 */
function phCytoscape (options) {
	// Create cytoscape instance in usual way
	let cy = cytoscape(options);
	// Method to (re)set cy.uobIdCounter
	cy.phSetIdCounter =  phSetIdCounter;
	// Counter used by functions in this document for ids when adding to cy
	cy.phIdCounter = cy.phSetIdCounter();
	// Other methods associated with our cytoscape instance (and not expected to be called from elsewhere)
	cy.phNextId = phNextId;
	cy.phListFunctions = phListFunctions;

	// Methods from P. Fowler-Wright
	// Section A - Explicit binding to cytoscape instance unnecessary
	cy.phHasGraphRealisation = phHasGraphRealisation;			// No binding necessary as doesn't expected a cytoscape instance
	cy.phHasBipartiteRealisation = phHasBipartiteRealisation; //  " -- "
	cy.phHasTreeRealisation = phHasTreeRealisation;
  cy.phTreesIsomorphic = phTreesIsomorphic;
	// Section B - both context (i.e. 'this') AND first argument of function is set as cy
	// This means each function can still be used in the normal way (e.g. myFunc(cy))
	 // First argument of bind() sets this; second fixes first argument of function. See https://javascript.info/bind
	cy.phMakeGraph = phMakeGraph.bind(cy, cy);
	cy.phMakeBipartite = phMakeBipartite.bind(cy, cy);
	cy.phMakeTree = phMakeTree.bind(cy, cy);
  cy.phMakeComplement = phMakeComplement.bind(cy, cy);
  cy.phDegreeSequence = phDegreeSequence.bind(cy, cy);
	cy.phEnumerateTrees = phEnumerateTrees.bind(cy, cy);
	return cy;
}

// General methods for our pseudo class cytoscape object 
/*
 * Increment & return the id counter used by functions in this document (to avoid id conflicts)
 * @returns {integer} incremented value of phIdCounter (bound to this === cy)
 */
function phNextId () {
	return ++this.phIdCounter;
}

/*
 * Adjust phIdCounter to avoid id conflicts with any elements already in cy (e.g. created via the 'options' parameter)
 * @returns {integer} maxID - integer to be assigned to cy.phIdCounter
 */
function phSetIdCounter () {
	// If integers have been used for ids already, return value greater than the maximum of these
	let maxID = -1;
	let idToNumber;
	// this refers to cy
	this.elements().forEach( ele => {
		idToNumber = Number(ele.id());
		maxID = Number.isSafeInteger(idToNumber) && idToNumber > maxID ? idToNumber : maxID;
	})
	return ++maxID; // 0 if no integer ids already in use
}

/*
 * Return array of (public) methods names bound to this
 */
function phListFunctions () {
	return Object.getOwnPropertyNames(this).filter(function (p) {
	    return typeof this[p] === 'function';
	});
}

// Functions from P. Fowler-Wright
// Section A - Not necessarily called via a phCytoscape object
/*
 * Tests whether a given array is consists of non-negative integers only. Used by many of my functions
 * @param {array} sequence - Array of values to be tested
 * @returns {bool} True if array consists of non-negative integers only.
 */
function phNonNegativeSequence (sequence) {
	// Check all elements safe integers (e.g. NOT Infinity, NaN) and positive
	// findIndex returns -1 if test function returns false for each element
	// So findIndex( ele => !Number.isSafeInteger(ele) || ele < 0 ) returns -1 if each element is a safe int and > 0
	return (sequence.findIndex( ele => !Number.isSafeInteger(ele) || ele < 0 ) == -1);
}
/*
 * Tests whether a sequence of non-negative integers is graphic using the Erdos-Gallai theorem
 * @param {array} degreeSeq - Array of vertex degrees. Needn't be ordered.
 * @returns {bool} Whether sequence is graphic or not
 */
function phHasGraphRealisation (degreeSeq) {
	// Check sequence contains non-negative integers only
  	if (!phNonNegativeSequence(degreeSeq)) { 
  		console.error('Vertex degrees must be non-negative.');
  		return false;
  	}
  	// Sort in reverse order as per theorem statement
  	degreeSeq = degreeSeq.sort((a,b) => b - a);
  	// Loop to calculate theorem formula at each k = 1...n where n is number of vertices
  	for (let k = 1; k <= degreeSeq.length; k++) {
  		// Get first k (highest degree) terms
  		let kHighestDegrees = degreeSeq.slice(0, k);
  		// And remaining terms
  		let nMinusKLowestDegrees = degreeSeq.slice(k);
  		// Sum over all k highest degrees = LHS of formula
  		let degreeSum = kHighestDegrees.reduce((sum, degree) => sum + degree, 0);
  		// Sum over min(d_i,k) for remaining degrees
  		let minSum = nMinusKLowestDegrees.reduce((sum, degree) => sum + Math.min(degree, k), 0);
  		// if, at any k, degreeSum > k(k-1 ) + minSum, sequence is non-graphic
  		if (degreeSum > k*(k-1) + minSum) {
  			return false;
  		}
  	}
  	// At all k, degreeSum <= k(k-1) + minSum <=> sequence is graphic
  	return true;
}
/*
 * Tests whether two sequences of non-negative integers have a simple bipartite graph realisation (Gale–Ryser theorem)
 * @param {array} degreeSeq1 - Integers describing the degrees of the first subgraph. Needn't be ordered.
 * @param {array} degreeSeq2 - Integers describing the degrees of the second subgraph Needn't be ordered.
 * @returns {bool}
 */
function phHasBipartiteRealisation (degreeSeq1, degreeSeq2) {
	// Check sequences contains positive integers only
  	if (!phNonNegativeSequence(degreeSeq1) || !phNonNegativeSequence(degreeSeq2) || degreeSeq1.includes(0) || degreeSeq2.includes(0)) { 
  		console.error('Vertex degrees of a bipartite graph must be positive.');
  		return false;
  	}
  	// Check sums of degrees are equal
  	if ( degreeSeq1.reduce((sum, degree) => sum + degree, 0) !== degreeSeq2.reduce((sum, degree) => sum + degree, 0)) {
  		console.log('Sums of degrees are unequal.');
  		return false;
  	}
  	// Use copies of degree sequences so as to avoid side-effects
  	let degreeSeq1Copy = degreeSeq1.slice();
  	let degreeSeq2Copy = degreeSeq2.slice();
  	// 'Stronger' version of Gale-Ryer thm. proved by Bender 2013 (see https://en.wikipedia.org/wiki/Gale%E2%80%93Ryser_theorem)
  	// To apply must have degree sequences of the same length and one (say the first) must be non-degreasing (i.e. 'reverse ordered')
  	let length1 = degreeSeq1Copy.length;
  	let length2 = degreeSeq2Copy.length;
  	if (length1 < length2) {
  		for (let i=0; i < length2 - length1; i++) {
	  		degreeSeq1Copy.push(0);
  		}
  	} else if (length1 > length2) {
  		for (let i=0; i < length1 - length2; i++) {
	  		degreeSeq2Copy.push(0);
  		}
  	}
  	degreeSeq1Copy = degreeSeq1Copy.sort( (a,b) => b - a );
  	// (Updated) length of each degree sequence, 'n'
  	let length = degreeSeq1Copy.length;
  	// The theorem: bigraphic iff an inequality holds for certain k=0...n-2 (case k = n-1 captured by check for equal sums)
  	for (let k=0; k < length - 1; k++) {
  		// No check necessary unless a_k > a_{k+1}
  		if ( degreeSeq1Copy[k] <= degreeSeq1Copy[k+1]) {
  			continue; 
  		}
  		// One side of inequality: sum of first (k+1) terms of degreeSeq1Copy
  		let degreeSeq1CopyTruncated = degreeSeq1Copy.slice(0,k+1);
  		let lhs = degreeSeq1CopyTruncated.reduce( (sum, degree) => sum + degree, 0);
  		// Other side: sum over second sequence of min(b_i, k+1) where b_i are terms of degreeSeq2Copy
  		let rhs = degreeSeq2Copy.reduce( (sum, degree) => sum + Math.min(degree, k+1), 0);
  		// Must have lhs <= rhs 
  		if ( lhs > rhs ) {
  			console.log('Sequences failed test of the Gale Ryser theorem (at k = ' + k + ').');
  			return false;
  		}
  	}
  	return true;
}
/*
 * Tests whether a sequence of natural numbers has a tree realisation
 * @param {array} degreeSeq - Array of vertex degrees. Needn't be ordered.
 * @returns {bool}
 */
function phHasTreeRealisation (degreeSeq) {
	// Check sequence contains non-negative integers only (and no invalid characters)
  	if (!phNonNegativeSequence(degreeSeq) || degreeSeq.includes(0)) { 
  		console.error('Vertex degrees of a tree must be non-negative.');
  		return false;
  	}
  	// The test: Sum of positive terms equals 2(n-2) <=> Has tree realisation (n is number of vertices)
  	return degreeSeq.reduce( (sum, degree) => sum + degree, 0) === 2*(degreeSeq.length - 1) ? true : false; 
}

// Section B - Expected to be called via a phCytoscape object

// Classes to conveniently encapsulate vertices and edges of a graph - used by phMakeGraph, phMakeBipartite, phMakeTree, phEnumerateTrees
// Each vertex has a unique id and a degree
class pVertex {
	constructor(degree) {
		this.id = cy.phNextId();
		this.degree = degree;
	}
}
// Each edge has an id as we as a 'source' and a 'target' (ids of the two vertices incident to the edge)
class pEdge {
	constructor(source, target) {
		this.id = cy.phNextId();
		this.source = source;
		this.target = target;
	}
}
/*
 * Calculate the degree sequence of a graph
 * @param {cytoscape} cytoscape - collection or instance containing cytoscape nodes
 * @param {bool} ascending -  whether the degree sequence should be in ascending (true) or descending (false) order
 * @returns {array} Ordered degree sequence 
 */
function phDegreeSequence (cy, ascending=true) {
	let degreeSeq = [];
	// Get nodes from degree (could be entire cytoscape instance)
	let nodes = cy.nodes();
	// Add degree of each node to list
	nodes.forEach( node => degreeSeq.push(node.degree()));
	// Return list, sorted according to whether ascending is true or not
	let sortSign = ascending ? 1 : -1;
	return degreeSeq.sort( (a, b) => sortSign*(a - b) );
}
/*
 * Perform Havel-Hakimi algorithm to generate a graph from a degree sequence
 * @param {cytoscape} cy - cytoscape instance to add graph to
 * @param {array} degreeSeq - Array of vertex degrees. Needn't be ordered
 * @param {bool} degreeLabels - Whether to label nodes with their degree
 * @param {bool} fixed - if true graph built will always be the same. N.B. The fixed procedure is such that the 
 constructed graph will be connected, IF a connected realisation exists.
 * @returns {cytoscape collection} - collection of nodes & edges added to cy
 */
function phMakeGraph (cy, degreeSeq, degreeLabels=false, fixed=true) {
	// Check sequence contains non-negative integers only
  	if (!phNonNegativeSequence(degreeSeq)) {
  		console.error('Vertex degrees must be non-negative.');
  		return false;
  	}
  	// If sum of degrees exceeds |E(K_n)| = n(n-1)/2 (total number of vertices in complete graph of n vertices),
  	// more efficient to perform algorithm on complement
  	let take_complement = false;
  	let numVerts = degreeSeq.length;
  	let double_num_edges = degreeSeq.reduce((sum, degree) => sum + degree, 0);
  	let number_complete = (numVerts)*(numVerts - 1)/2;
  	if (double_num_edges > number_complete) {
  		// Degree d_i maps to (n-1) - d_i under complement
  		degreeSeq = degreeSeq.map( degree => (numVerts - 1) - degree);
  		// After performing the algorithm we will need to take the complement
  		take_complement = true;
  	}
 	// Flag to indicate whether sequence is graphic (true) or not (false); set on termination of algorithm
	let graphic;
	// Array of vertices describing graph
	let vertArr = [];
	degreeSeq = degreeSeq.sort( (a,b) => b - a ); // Decreasing order
	degreeSeq.forEach( degree => vertArr.push(new pVertex(degree)) );
  	// Array of edges used to create connected graph in cytoscape, if cy provided. 
  	let edgeArr = []; 
  	// Perform the recursive algorithm on a copy of vertArr; sets graphic and constructs edgeArr
	let vertArrCopy = JSON.parse(JSON.stringify(vertArr));
  	havelHakimi(vertArrCopy);
  	// If cy == null, just return whether sequence is graphic
  	if (cy == null) {
  		return graphic;
  	}
  	// Otherwise we construct a graph (if possible) using edgeArr and add to cytoscape instance
  	if (!graphic) {
  		console.error('Degree sequence is non-graphic.');
  		return false;
  	}
  	// If the 'complement degrees' were used above now must take complement of graph represented by edgeArr
  	if (take_complement) {
  		let complementEdgeArr = [];
  		for (let i=0; i < numVerts; i++) {
  			for (let j=0; j < i; j++) {
  				if (edgeArr.findIndex( edge => 
  					edge.source === vertArr[i].id && edge.target === vertArr[j].id ||
  					edge.source === vertArr[j].id && edge.target === vertArr[i].id
  					) === -1) {
  					complementEdgeArr.push(new pEdge(vertArr[i].id, vertArr[j].id));
  				}
  			}
  		}
  		edgeArr = complementEdgeArr;
  	}
  	// Collection of elements added by this function, returned to caller 
  	// Also allows us to edit the style/layout for added graph without affecting rest of cytoscape instance.
  	let havelCollection = cy.collection();
  	buildGraph();
	// Pan and zoom to fit all elements of cytoscape instance
  	cy.fit();
	// cy.fit(havelCollection); // Pan and zoom to newly created graph (only)
  	return havelCollection;


  	// ~ Description of procedure (fixed=true) ~
  	/* Choose random (final) vertex sequence and 'attach' to the d vertices remaining vertices with the largest degrees.
  	Here d is the degree of the chosen vertex and 'attach' translates to: Create an pEdge, add it to edgeArr, and
  	decrement target vertex degrees. Remove the chosen vertex as well as any other vertex whose degree is turned 0. 
  	Repeat until all vertices exhausted or algorithm fails due to inability to 'exhaust' chosen node.
  	This is performed recursively.
  	*/
  	function havelHakimi (remainingVerts) {
  		// Remove any exhausted vertices (done prior to .pop() as if take_complement may have vertices of 0 degree)
  		remainingVerts = remainingVerts.filter( vertex => vertex.degree !== 0);
  		// ALgorithm continues until all vertices exhausted (positive solution)
  		if (remainingVerts.length === 0) {
  		  	graphic = true;
  		  	return;
  		}
  		// Sort the remanding array of vertices in descending order of degree
  		remainingVerts = remainingVerts.sort( (v1, v2) => v2.degree - v1.degree );
  		// If want to generate same graph each time, always select vertex with smallest degree. Otherwise choose at random
		let chosenVertex;
		if (fixed) {
			// N.B. Choosing smallest degree vertex guarantees a connected graph, if such a realisation exists
			// See "Building connected graphs" http://szhorvat.net/pelican/hh-connected-graphs.html
			chosenVertex = remainingVerts.pop();
		} else {
			chosenVertex = remainingVerts.splice(Math.floor(Math.random()*remainingVerts.length), 1)[0];
		}
		// One iteration: Remove vertex of (smallest) degree X; subtract 1 from the degrees of vertices with X largest degrees
		// e.g. [5,4,4,2] -> [4,3,4]
		// Algorithm fails (graph non graphic) if not enough vertices to connect to
		if (chosenVertex.degree > remainingVerts.length) {
			graphic = false;
			return;
		}
		for (let i=0; i < chosenVertex.degree; i++) {
			// Subtract 1 from degree number of remaining vertices, starting with that which has the largest degree
			remainingVerts[i].degree -= 1;
			// Add an edge between chosen and target vertex
			edgeArr.push(new pEdge(chosenVertex.id, remainingVerts[i].id)); 
		}
		// Continue until all vertices exhausted
		return havelHakimi(remainingVerts);
  	}

  	function buildGraph () {
		// Add vertices (cytoscape nodes)
		vertArr.forEach ( vertex => {
			let node = cy.add({ group: 'nodes', data: { id: vertex.id }});
			havelCollection = havelCollection.union(node);
		});
		// Add edges generated by addLeaf recursive procedure. 
		edgeArr.forEach( edge => {
			let cyEdge = cy.add([{ group: 'edges', data: { id: edge.id, source: edge.source, target: edge.target } }]);
			havelCollection = havelCollection.union(cyEdge);
		});
		// Label each vertex according to its degree, if requested
		if (degreeLabels) {
			let nodeCollection = havelCollection.nodes();
			nodeCollection.forEach( (node) => node.style('label', (ele) => ele.degree()) );
		}
		// Set and run a layout for the graph  (affects havelCollection only)
		var layout = havelCollection.layout({
			name: 'cose', // Compound Spring Embedder layout -  uses a physics simulation to lay out graphs.
			animate: false // Don't animate when running the layout
		});
		// run layout
		layout.run();
  	}
}
/*
 * Construct a bipartite graph from two degree sequences, if possible
 * @param {cytoscape} cy - cytoscape instance to add graph to
 * @param {array} degreeSeq1 - Array of degrees for first set of vertices
 * @param {array} degreeSeq2 - Array of degrees for second set of vertices
 * @param {bool} degreeLabels - Whether to label nodes with their degree
 * @param {bool} fixed - if true graph built will always be the same
 * @returns {array} - Three cytoscape collections containing, respectively the nodes of the 1st, the 2nd subgraph, and the edges they share
 */
function phMakeBipartite (cy, degreeSeq1, degreeSeq2, degreeLabels=false, fixed=false) {
	if (!phHasBipartiteRealisation(degreeSeq1,degreeSeq2)) {
		console.error('Degree sequences are not bipartite.');
		return false;
	}
	// Order sequences (means order of sequences passed to this function does not matter when fixed=true)
	degreeSeq1 = degreeSeq1.sort( (a,b) => b - a );
	degreeSeq2 = degreeSeq2.sort( (a,b) => b - a );
	// Array of vertex objects for each subgraph
	vertArr1 = [];
	vertArr2 = [];
	degreeSeq1.forEach( degree => vertArr1.push(new pVertex(degree)) );
	degreeSeq2.forEach( degree => vertArr2.push(new pVertex(degree)) );
	// Array of edges between subgraphs
	edgeArr = [];
	// Initialisation complete, create graphs using deep copies of vertArr1, vertArr2
	let vertArr1Copy = JSON.parse(JSON.stringify(vertArr1));
	let vertArr2Copy = JSON.parse(JSON.stringify(vertArr2));
	addVertex(vertArr1Copy, vertArr2Copy);
	// Construct 3 cytoscape collections describing subgraphs and edges between them
	let bipartiteCollection1 = cy.collection();
	let bipartiteCollection2 = cy.collection();
	let edgeCollection = cy.collection();
	let entireCollection = cy.collection();
	buildGraph();
	// Pan and zoom to fit all elements of cytoscape instance
	cy.fit();
	// cy.fit(entireCollection); // Same but for constructed graph only
	//  return graph as array of cytoscape collections: [graph1, graph2, edges]
	return [bipartiteCollection1,bipartiteCollection2,edgeCollection];

	// ~ Description of procedure (fixed=true) ~
	/* Choose random (first) vertex from first sequence and 'attach' to the d vertices of the second sequence with the
	LARGEST degree. Here d is the degree of the chosen vertex and 'attach' translates to: Create an pEdge, add it to 
	edgeArr, and decrement target vertex degrees. Remove the chosen vertex as well as any vertex in second sequence
	whose degree is turned 0. Repeat until all vertices exhausted. This is performed recursively.
	*/
	function addVertex(remainingVerts1, remainingVerts2) {
		// Process stops when all vertices exhausted (from both graphs)
		if (remainingVerts1.length === 0) {
			return;
		}
		// Sort remaining vertices of second sequence according to degree (decreasing)
		remainingVerts2 = remainingVerts2.sort( (v1, v2) => v2.degree - v1.degree );
		// Choose vertex randomly unless fixed=true (for greater variability could also choose target randomly if multiple vertices with same degree)
		let chosenSourceIndex = 0;
		let chosenTargetIndex = 0;
		if (!fixed) {
			chosenSourceIndex = Math.floor(Math.random() * remainingVerts1.length);
			// Target must have largest degree from remaining vertices (remainingVerts2[0] is always of maximum degree due to ordering)
			let numberMaximumDegreedVertices = remainingVerts2.findIndex( vertex => vertex.degree !== remainingVerts2[0].id); 
			chosenTargetIndex = Math.floor(Math.random() * numberMaximumDegreedVertices);
		}
		// Remove chosen source from remainingVerts1
		let source = remainingVerts1.splice(chosenSourceIndex,1)[0]; // N.B. splice returns an ARRAY - select first (only) element
		// Decrement degrees of targets whilst adding edges to edgeArr
		for (let i=0; i < source.degree; i++) {
			// Aside: For enumeration problem would have to iterate through all possible choices of targets, e.g.
			// if source.degree=1 and second sequence is (3,3,3,2,2,1), i = 0, 1, 2 are all valid choices
			remainingVerts2[i].degree -= 1;
			edgeArr.push(new pEdge(source.id, remainingVerts2[i].id));
		}
		// Remove any exhausted vertices from second array
		remainingVerts2 = remainingVerts2.filter( vertex => vertex.degree !== 0);
		// Continue until all vertices exhausted
		addVertex(remainingVerts1, remainingVerts2);
	}

	function buildGraph() {
		// Add vertices (cytoscape 'nodes')
		vertArr1.forEach ( vertex => {
			let node = cy.add({ group: 'nodes', data: { id: vertex.id }});
			bipartiteCollection1 = bipartiteCollection1.union(node);
		});
		vertArr2.forEach ( vertex => {
			let node = cy.add({ group: 'nodes', data: { id: vertex.id }});
			bipartiteCollection2 = bipartiteCollection2.union(node);
		});
		// Add edges generated by addLeaf recursive procedure. 
		edgeArr.forEach( edge => {
			let cyEdge = cy.add([{ group: 'edges', data: { id: edge.id , source: edge.source, target: edge.target } }]);
			edgeCollection = edgeCollection.union(cyEdge);
		});
		entireCollection = bipartiteCollection1.union(bipartiteCollection2).union(edgeCollection);
		// Label each vertex according to its degree, if requested
		if (degreeLabels) {
			// All nodes
			let nodeCollection = entireCollection.nodes();
			nodeCollection.forEach( (node) => node.style('label', (ele) => ele.degree()) );
		}
		// Set and run a layout for the graph
		var layout = entireCollection.layout({
			name: 'cose', // Compound Spring Embedder layout -  uses a physics simulation to lay out graphs.
			animate: false // Don't animate when running the layout
		});
		layout.run();
	}
}

/*
 * Construct a tree from a given degree sequence, if possible
 * @param {cytoscape} cy - cytoscape instance to add graph to
 * @param {array} degreeSeq - Array of vertex degrees. Needn't be ordered.
 * @param {bool} degreeLabels - Whether to label nodes with their degree
 * @param {bool} fixed - if true tree built will always be the same
 * @returns {cytoscape collection} - collection of nodes & edges added to cy
 */
function phMakeTree (cy, degreeSeq, degreeLabels=false, fixed=false) {
	// Check sequence has a tree realisation
	if (!phHasTreeRealisation(degreeSeq)) {
		console.error('Degree sequence has no tree realisation.');
		return false;
	}
	// Construct array of Vertex objects for the tree. For simplicity firstly reverse order degreeSeq (largest -> smallest)
	// so ids 0,1,... are assigned in order of decreasing vertex degree. A distinct copy is made which is mutated in the below procedure.
	let vertArr = [];
	degreeSeq = degreeSeq.sort( (a,b) => b - a ); // Decreasing order
	degreeSeq.forEach( degree => vertArr.push(new pVertex(degree)) );
	// Make a deep copy (JSON serialization method works provided no functions or complex types in objects)
	let vertArrCopy = JSON.parse(JSON.stringify(vertArr));
	// Array of Edge objects for tree, constructed in below while loop
	let edgeArr = [];
	// Initialisation complete, initiate recursive procedure to populate edgeArr
	addLeaf(vertArrCopy);
	// Construct cytoscape collection describing tree
	let treeCollection = cy.collection();
	buildTree();
	// Pan and zoom to fit all elements of cytoscape instance
	cy.fit(); 
	// ct.fit(treeCollection) // Centre on constructed graph only
	return treeCollection;


	// ~ Description of Procedure ~
	/* 
	Take any leaf (from vertArr) and attach to a branch (make an pEdge object), removing the leaf
	from vertArr and subtracting 1 from the degree of the chosen branch. Repeat with the modified vertArr (achieved with
	a recursive call), adding the pEdge at each step to a growing array called edgeArr (initially empty []). 
	Stop when the degree sequence remaining is just (1,1) (length 2) - the final edge is added between the two corresponding vertices.
	*/
	function addLeaf(remainingVerts) {
		if (remainingVerts.length == 2) {
			// Final edge between those last two 'leaves'
			edgeArr.push(new pEdge(remainingVerts[0].id, remainingVerts[1].id)); // Recall pEdge has a 'source' and a 'target' 
			return;
		}
		// Sort remaining vertices according to degree (decreasing)
		remainingVerts = remainingVerts.sort( (v1, v2) => v2.degree - v1.degree );
		// Choose leaf and branch randomly unless fixed=true in which case take final leaf and first branch (greatest degree)
		let chosenLeafIndex;
		let targetBranchIndex;
		if (fixed) {
			chosenLeafIndex = remainingVerts.length-1;
			targetBranchIndex = 0;
		} else {
			// Index of first leaf. Always > 2 as final case of (1,1) caught by above if () condition.
			let firstLeafIndex = remainingVerts.findIndex( vertex => vertex.degree === 1 );
			let numberLeaves = remainingVerts.length - firstLeafIndex;
			let numberBranches = firstLeafIndex;
			chosenLeafIndex = firstLeafIndex + Math.floor(Math.random() * numberLeaves);
			targetBranchIndex = Math.floor(Math.random() * numberBranches);
		}
		// Remove chosen leaf from remainingVerts
		let chosenLeaf = remainingVerts.splice(chosenLeafIndex,1)[0];
		// Decrement degree of chosen branch
		remainingVerts[targetBranchIndex].degree -= 1;
		// Add edge between chosen leaf and branch
		edgeArr.push(new pEdge(chosenLeaf.id, remainingVerts[targetBranchIndex].id));
		// Continue until all vertices exhausted
		addLeaf(remainingVerts);
	}

	function buildTree() {
		// Add vertices (cytoscape 'nodes')
		vertArr.forEach ( vertex => {
			let node = cy.add({ group: 'nodes', data: { id: vertex.id }});
			treeCollection = treeCollection.union(node);
		});
		// Add edges generated by addLeaf recursive procedure. Use 'e' + edge.source + '-' +  edge.target as id
		edgeArr.forEach( edge => {
			let cyEdge = cy.add([{ group: 'edges', data: { id: 'e' + edge.source + '-' + edge.target , source: edge.source, target: edge.target } }]);
			treeCollection = treeCollection.union(cyEdge);
		});
		// Label each vertex according to its degree, if requested
		if (degreeLabels) {
			let nodeCollection = treeCollection.nodes();
			nodeCollection.forEach( (node) => node.style('label', (ele) => ele.degree()) );
		}
		// Set and run a layout for the graph
		var layout = treeCollection.layout({
			name: 'cose', // Compound Spring Embedder layout -  uses a physics simulation to lay out graphs.
			animate: false // Don't animate when running the layout
		});
		layout.run();
	}
}

/*
 * Attempts to enumerate all trees with a given degree sequence, up to isomorphism (i.e. only non-isomorphic trees are returned)
 * @param {array} degreeSeq - Array of vertex degrees. Needn't be ordered.
 * @returns {array} - First element is a cytoscape collection describing the nodes of any of the possible trees. Following
 * elements are collections of edges, each describing a (different) tree with vertices described. 
 * Nothing is added to cy; to draw graph use cy.add() on the node collection followed by any one of the edge collections
 * N.B. This is a brute force procedure with repeated calls to a cytoscape instance and so not efficient by any means
 * The algorithm used to check for isomorphism - AHU - is also not the most efficient). It is also not well tested.
 *
 * //Example code using the return of this function (assuming let cy = phCytoscape(...) has been initialised):
 * let collections = cy.phEnumerateTrees([1,1,1,1,3,3]); // Simple example sequence; 6 isomorphic trees are generated by this function
 * cy.add(collections[0]); // Draw the nodes
 * console.log(collections.length - 1) // Check how many trees were generated. OUTPUT: 6
 * cy.add(collections[3]); // Draw 3rd set of edges (for example)
 * // Now make an run a layout to display graph nicely
 * let layout = cy.layout({
 * 	name: 'cose',
 *	animate: false
 *	});
 * layout.run();
 */
function phEnumerateTrees(cy, degreeSeq) {
	// Check sequence has a tree realisation
	if (!phHasTreeRealisation(degreeSeq)) {
		console.error('Degree sequence has no tree realisation.');
		return false;
	}
	// Holds edges of all generated trees for the given degree sequence (Array of arrays, each array the edges of one possible tree)
	let allEdges = [];
	// Construct array of pVertex objects for the tree. For simplicity firstly reverse order degreeSeq (largest -> smallest)
	// so ids 0,1,... are assigned in order of decreasing vertex degree. 
	let vertArr = [];
	degreeSeq = degreeSeq.sort( (a,b) => b - a ); // Decreasing order
	degreeSeq.forEach( degree => vertArr.push(new pVertex(degree)) );
	// Initialisation complete, now populate allEdges using a recursive procedure
	attachLeaf(vertArr);
	// Remove duplicates from allEdges. Note this does NOT check for isomorphism
	removeDuplicates();
	// Array of cytoscape collections to be returned; first element a collection of nodes, all others collection of edges
	let rtnArr = [];
  // Cytoscape collections of edges, to be pruned for isomorphism
  let cyEdgeArr = [];
  // Empty cytoscape object to hold the collections (Note no 'container' option hence no display of these elements)
  let tempCy = cytoscape();
  // Add collection of nodes to rtnArr and fill cyEdgeArr
	genCyCollections();
  // Removes duplicates up to isomorphism from cyEdgeArr
  removeIsomorphic();
  // Append remaining edge collections onto rtnArr
  rtnArr.push(...cyEdgeArr);
	// Notify how many sets of edges were generated
	console.log('Edges for ' + (rtnArr.length - 1) + ' trees generated (indices 1..' + (rtnArr.length - 1) + ' of returned array).');
	return rtnArr;

	// ~ Description of Procedure ~
	/* 
	To construct ONE tree: Take any leaf (from vertArr) and attach to a branch (make an pEdge object), removing the leaf
	from vertArr and subtracting 1 from the degree of the chosen branch. Repeat with the modified vertArr (achieved with
	a recursive call), adding the Edge at each step to a growing array called edgeArr (initially empty []). 
	Stop when the degree sequence remaining is just (1,1) - the final edge is added between the two corresponding vertices.
		
	To record ALL trees which may be constructed in this manner at EACH step we iterate through ALL the possible choices
	of leaf ('i' for-loop below) and all the possible choices of branch ('j' for-loop below). for each choice at every 
	step a separate call to attachLeaf (creating a distinct copy of the edge arrays) is made so the process continually 
	'branches' out, generating many different arrays of edges, each detailing a possible tree with the vertices in vertArr/
	*/ 
	function attachLeaf (remainingVerts, edgeArr=[]) {
		// As stated  above, the construction finishes when we have been reduced to an effective degree sequence of (1,1) (length 2)
		if (remainingVerts.length == 2) {
			// Create a new array of Edges which is a copy of the array from the last step (if we just added to the original array,
			// the edges from the trees of all possible choices would be added to the same array!).
			//  This is a 'deep' copy - important as elements of array of objects, not basic data types
			let newEdgeArray = JSON.parse(JSON.stringify(edgeArr));  // N.B. This JSON shortcut only works if no functions or complex types 
			// Final edge between those last two 'leaves'
			newEdgeArray.push(new pEdge(remainingVerts[0].id, remainingVerts[1].id)); // Recall pEdge has a 'source' and a 'target' 
			// Edges of tree complete so add to allEdges; time to head back up the recursive tree!
			allEdges.push(newEdgeArray); 
			return;
		}
		// Sort the remanding array of vertices in descending order of degree (makes selection of leaf & branch easy)
		remainingVerts = remainingVerts.sort( (v1, v2) => v2.degree - v1.degree );
		// Index of first leaf in remaingVertices; all higher indices are definitely leaves due to sorting
		// N.B. firstLeafIndex >= 1 (the final case of a degree sequence of (1,1) is caught in the above if ())
		let firstLeafIndex = remainingVerts.findIndex( vertex => vertex.degree === 1 );
		// Number of leaves and branches we have to iterate over to explore all possible constructions
		let numberLeaves = remainingVerts.length - firstLeafIndex;
		let numberBranches = firstLeafIndex;
		// Iterating over the possible leaf choices
		for (let i=0; i < numberLeaves; i++) {
			// Make a deep copy of remainingVerts (otherwise the same array would be mutated for successive values of i)
			let remainingVertsCopy = JSON.parse(JSON.stringify(remainingVerts));
			// Remove chosen leaf from remainingVertsCopy
			let chosenLeaf = remainingVertsCopy.splice(firstLeafIndex + i,1)[0];
			// Iterate over possible branch choices, given a leaf choice
			for (let j=0; j < numberBranches; j++) {
				// We need another deep copy (otherwise different values of j would modify the same objects)
				let remainingVertsCopy2 = JSON.parse(JSON.stringify(remainingVertsCopy));
				// Subtract 1 from the degree of the chosen branch
				let targetBranch = remainingVertsCopy2[j];
				targetBranch.degree -= 1;
				// Create deep copy of edgeArr (each 'branch' of this process should result in its own array of edges)
				let newEdgeArray = JSON.parse(JSON.stringify(edgeArr));
				// Add edge between chosen leaf and branch, then call attachLeaf with the modified array of vertices
				newEdgeArray.push(new pEdge(targetBranch.id, chosenLeaf.id));
				attachLeaf(remainingVertsCopy2, newEdgeArray);
			}
		}
	}

	function removeDuplicates () {
		// Order each array of edges to make comparisons between arrays easier
		// 1) order the ids of each individual edge - ensure 'source' id is always less than or equal to 'target id
		allEdges.forEach( edgeArr => {
			let source;
			edgeArr.forEach ( edge => {
				// Save id of source
				source = edge.source;
				if (source > edge.target) {
					// Swap source and target ids
					edge.source = edge.target;
					edge.target = source;
				}
			})
		});
		// 2) order the edges according to the source id, from smallest to largest source id
		allEdges.forEach( edgeArr => edgeArr.sort( (edge1, edge2) => edge1.source - edge2.source || edge1.target - edge2.target) );
		/* Example: edgeArr = [ [1,0], [5,2], [331] ] ( ['source', 'target'])
		1) -> [ [0,1], [2,5], [1,3] ] 
		2) -> [ [0,1], [1,3], [2,5] ] 
		*/
		// Array of indices (of allEdges) flagging duplicates
		let duplicateIndicies = [];
		// Procedure: Take 1st element of allEdges. Compare to 2nd, 3rd,.. entries, adding to duplicateIndicies is a match is found
		// Next take 2nd element of allEdges and compare to 3rd, 4th,...
		// At any point if an edge has already been added to the duplicate list, skip it (continue's below).
		for (let i=0; i < allEdges.length; i++) {
			if (duplicateIndicies.indexOf(i) != -1) {
				// allEdges[i] has already been marked as a duplicate
				continue;
			}
			// N.B. j starts at i+1
			for (let j=i+1; j < allEdges.length; j++) {
				if (duplicateIndicies.indexOf(j) != -1) {
					// allEdges[j] has already been marked as a duplicate
					continue;
				}
				// Compare the two arrays of edges; if identical add j to duplicateIndices
				if (treesIdentical(allEdges[i], allEdges[j])) {
					duplicateIndicies.push(j);
				}
			}
		}
		// Finally, remove all duplicates
		allEdges = allEdges.filter( (tree, index) => duplicateIndicies.indexOf(index) === -1 );
	}

	function treesIdentical(edgeArr1, edgeArr2) {
		// To test whether two arrays of edges are identical, can test the source and edge ids of successive terms
		// This works because the arrays have been ordered as described in removeDuplicates
		for (let i=0; i < edgeArr1.length; i++) {
			if (edgeArr1[i].source != edgeArr2[i].source || edgeArr1[i].target != edgeArr2[i].target) {
				// There is a difference
				return false;
			}
		}
		// No difference found between two arrays
		return true;
	}

	// Generate cytoscape collections from the arrays of pVertex and pEdge arrays
	function genCyCollections () {
		// Collection of nodes; these can go directly onto rtnArry
	 	let nodeCollection = tempCy.collection();
		vertArr.forEach ( vertex => {
			let node = tempCy.add({ group: 'nodes', data: { id: vertex.id }});
			nodeCollection = nodeCollection.union(node);
		});
		rtnArr.push(nodeCollection);
		allEdges.forEach( edgeArr => {
			// Collection of edges for one tree. These populate entries (index) 1,2,... or rtnArry
			let edgeCollection = tempCy.collection();
			edgeArr.forEach( edge => {
				let cyEdge = tempCy.add({ group: 'edges', data: { id: edge.id, source: edge.source, target: edge.target }});
				edgeCollection = edgeCollection.union(cyEdge);
			});
			cyEdgeArr.push(edgeCollection);
		} );
	}

  // Remove edges from cyEdgeArr which correspond to (given the fixed set of vertices) isomorphic trees
  // This is done via calls to phTreesIsomorphic
  function removeIsomorphic () {
    // Iterate backwards through cyEdgeArr: Start with final element and compare to first, second, third,...(final-1)
    // If a match (isomorphic tree) is found, remove the final element and continue with the reduced array
    for (let i=(cyEdgeArr.length - 1); i >= 0; i--) {
      // Obviously do not test a tree against itself (i===j)
      for (let j=0; j < i; j++) {
        // Cytoscape collections describing the trees to be tested (phTreesIsomorphic requires full tree, not just edges)
        let tree1 = rtnArr[0].union(cyEdgeArr[i]); // Recall rtnArr[0] is collection of nodes
        let tree2 = rtnArr[0].union(cyEdgeArr[j]);
        if (phTreesIsomorphic(tree1, tree2)) {
          // console.log('Edges at ' + j + ' & ' + i + ' describe isomorphic trees, deleting ' + i);
          // Remove duplicate. As iterating through array backwards this does not cause issues
          cyEdgeArr.splice(i,1); // Modifies in place
          // Break from inner loop
          break;
        }
      }
    }
  }
}

/*
 * Build the complement of a graph associated with a cytoscape instance 
 * @param {cytoscape instance} cy - cytoscape instance containing a graph
 * @param {bool} newGraph - If true make a new graph instead of using current nodes and removing previous edges. Nodes
 * of new graph are coloured a nice shade blue by default
 * @param {string} newGraphNodeColour - HTML Colour name/Hex Code used to colour complement nodes, if newGraph
 * @returns {cytoscape collection} - collection of nodes & edges describing complement of graph
 */
function phMakeComplement (cy, newGraph=false, newGraphNodeColour="MEDIUMSLATEBLUE") {
	// Get nodes and edges of current edges
	let nodes = cy.nodes();
	let edges = cy.edges();
	// Make list of node ids in graph
	let nodeIds = [];
	nodes.forEach( node => nodeIds.push( node.id() ) );
	// Empty collection to hold edges and nodes added to cy; returned by this function
	let complementCollection = cy.collection();
	// If making a new graph for complement, add nodes with same ids as current graph but a 'copy' prefix
	if (newGraph) {
		let complementNodes = [];
		nodes.forEach( node => complementNodes.push( { group: 'nodes', data: { id: 'copy' + node.id() }} ));
		complementCollection = cy.add(complementNodes);
		// Define and apply a style class to the complement nodes to colour them
		cy.style().selector('.complementNode').style().style('background-color', newGraphNodeColour).update();
		complementCollection.addClass('complementNode');
	} else {
		// Otherwise use current nodes and remove current edges
		complementCollection = nodes;
		cy.remove('edge');
	}
	// Edges (plain objects) of complement
	let complementEdges = [];
	// Iterate through each pair of nodeIds and see if any edge in current graph is incident to those edges
	// If not (edgeFound = false), edge is added to array completmentEdges
	for (let i=0; i < nodeIds.length; i++) {
		for (let j=0; j < i; j++) { // N.B. j < i (otherwise would see each id pairs twice, as well as j=i which would check for loops)
			let edgeFound = false;
			edges.forEach( edge => {	// N.B. Can't use findIndex method on cy.egdes() object.
				let sourceId = edge.source().id();
				let targetId = edge.target().id();
				// Check both possible orderings
				if (sourceId === nodeIds[i] && targetId === nodeIds[j] || sourceId === nodeIds[j] && targetId === nodeIds[i]) {
					edgeFound = true;
				}
				});
			if (!edgeFound) {
				// If building new graph, need to prefix source/target ids with 'copy' as above
				if (newGraph) {
					complementEdges.push({ group: 'edges', data: { id: cy.phNextId(), source: 'copy' + nodeIds[i], target: 'copy' + nodeIds[j] }});
				} else {
					complementEdges.push({ group: 'edges', data: { id: cy.phNextId(), source: nodeIds[i], target: nodeIds[j] }});
				}
			}
		}
	}
	// Add edges to the growing collection
	complementCollection = complementCollection.union(cy.add(complementEdges));
	// If making a new graph run a layout on entire graph (otherwise graphs may be placed on top one another etc.)
	if (newGraph) {
		let layout = cy.elements().layout({ // Alternatively could apply to complementCollection instead of cy.elements() 
			name: 'cose', // Compound Spring Embedder layout -  uses a physics simulation to lay out graphs.
			animate: false // Don't animate when running the layout
		});
		// run layout
		layout.run();
		// Fit to view all the nodes in the graph
		cy.fit(); 
	}
	return complementCollection;
}

/*
 * Test whether two (ordinary) trees are isomorphic
 * @param {cytoscape} cy1 - cytoscape instance (or collection) containing first tree
 * @param {cytoscape} cy2 - cytoscape instance (or collection) containing second tree
 * @returns {bool}
 */
function phTreesIsomorphic(cy1, cy2) {
  // Check each cytoscape instance describes a tree
  if (!phIsTree(cy1) || !phIsTree(cy2)) {
    console.error('One of graphs is not a tree.');
    return;
  }
  // console.log('Both graphs are trees.');
  // Find the centre(s) of each tree
  let t1Centre = phFindCentres(cy1);
  let t2Centre = phFindCentres(cy2);
  // Three cases
  // a) Trees have different number of centres -> not isomorphic
  if (t1Centre.length !== t2Centre.length) {
    // console.log('Trees have different count of centres hence are not isomorphic.');
    return false;
  }
  // b) each tree has only one centre
  if (t1Centre.length === 1) {
    // Perform rooted tree isomorphism check with centres as roots
    // console.log('Performing rooted isomorphism check with unique centres.');
    return phRootedTreesIsomorphic(cy1, t1Centre[0], cy2, t2Centre[0]);
  }
  // c) Each tree has two centres
  // Perform rooted tree check twice, changing the centre used as the root for the first tree in the second check
  // console.log('Performing rooted isomorphism check with two different centres:');
  return phRootedTreesIsomorphic(cy1, t1Centre[0], cy2, t2Centre[0]) || phRootedTreesIsomorphic(cy1, t1Centre[1], cy2, t2Centre[0]);
}
/*
 * Test whether a graph is a tree
 * @param {cytoscape} collection - cytoscape instance or collection containing the graph
 * @returns {bool}
 */
function phIsTree(collection) {
  let edges = collection.edges();
  let nodes = collection.nodes();
  // Check for loops
  let containsloops = false;
  edges.forEach( edge => {if (edge.isLoop()) {containsloops = true;}});
  if (containsloops) {
    // console.log('Graph contains loops (not simple).');
    return false;
  }
  // Check connected
  if (!phIsConnected(collection)) {
    // console.log('Not all vertices are connected.');
    return false;
  }
  // Any connected graph on n vertices with n-1 edges is a tree (R. Leek Algebra and combinatorics Lecture notes Corollary 5.22)
  if (edges.length === nodes.length - 1) {
    // console.log('Graph is connected and the number of edges equals number of vertices minus 1.');
    return true;
  } else {
    // console.log('Graph is connected but the number of edges does not equal number of vertices minus 1.');
    return false;
  }
}
/*
 * Test whether a graph is connected
 * @param {cytoscape} collection - cytoscape instance or collection containing the graph
 * @returns {bool}
 * The dijkstra routine from cytoscape.js is used to determine the shortest path from an arbitrarily chosen vertex (the first
 * one) and every other vertex in the graph; this returns Infinity if no path exists
 */
function phIsConnected(collection) {
  let nodes = collection.nodes();
  // .dijkstra returns an object with a 'distanceTo(node)' function which returns the path length from the root to node
  let dij = collection.dijkstra({root:'#'+nodes[0].id()}); // Use first node as root
  // Flag to indicated connected or not (cannot return from within a .forEach loop)
  let connected = true;
  nodes.forEach ( node => {
    if (dij.distanceTo('#' + node.id()) === Infinity) {
      // No path exists, so graph is not connected
      connected = false;
    }
  });
  return connected;
}
/*
 * Find the centre(s) of a tree
 * @param {cytoscape} collection - cytoscape instance or collection containing the tree
 * @returns {cytoscape} - collection containing centre or centres of tree
 * This is done iteratively: At each step remove all leaves from the tree, until only one or two vertices remain (the centre(s))
 */
function phFindCentres(collection) {
  // So as to not affect original graph we make a copy
  let cyCopy = cytoscape(); //  No container hence instance is not displayed
  cyCopy.add(collection.clone()); // clone collection and add to our hidden instance
  let remainingNodes = cyCopy.nodes(); // collection of remaining nodes
  let remaingNodeCount = remainingNodes.length; // Procedure ends when only one or two nodes left
  let toRemove = []; // Array of nodes to be removed at the end of this iteration
  while (remainingNodes.length > 2) {
    for (let i = 0; i < remainingNodes.length; i++ ) {
      // Mark all leaves for removal from the tree
      if (remainingNodes[i].degree() === 1) {
        // console.log('Node ' + remainingNodes[i].id() + ' is now a leaf; removing.');
        toRemove.push(remainingNodes[i]);
      }
    }
    // Remove all leaves (done AFTER the for loop). Note that this removes their edges from cyCopy too
    toRemove.forEach( node => cyCopy.remove(node) );
    // Clear toRemove and get the remaining nodes for the next iteration
    toRemove = [];
    remainingNodes = cyCopy.nodes();
  }
  return remainingNodes;
}

/*
 * Test whether two rooted trees are isomorphic using the AHU algorithm
 * @param {cytoscape} tree1 - cytoscape instance or collection containing the first tree
 * @param {cytoscape} root1 - cytoscape element (in tree1) that is the root of tree 1
 * @param {cytoscape} tree2 - cytoscape instance or collection containing the second tree
 * @param {cytoscape} root2 - cytoscape element (in tree2) that is the root of tree 2
 * @returns {bool}
 * The algorithm is described p.84 of Alfred V. Aho, John E. Hopcroft, Jeffrey D. Ullman | The Design and Analysis of Computer Algorithms
 */
function phRootedTreesIsomorphic (tree1, root1, tree2, root2) {
  // Class to encapsulate vertex information used in the algorithm 
  // Whilst cytoscape implements depth- and breadth-first searches, we work from level 0 up the tree (so don't directly use these search tools)
  class Vertex {
    constructor(id, depth, parentId) {
      // Unique Cytoscape idea
      this.id = id;
      // Depth of vertex is distance (length of path) from root
      this.depth = depth;
      // Cytoscape id of immediate parent
      this.parentId = parentId;
      // Tuple assigned to each branch (non-leaf) on a level
      this.assignedTuple = [];
      // Integer assigned to all vertices on a level
      this.assignedInt;
      // Height of tree is maximum over vertex depths
      this.treeHeight;
    }
    get level() {
      // Level is height of tree minus depth of his vertex 
      return this.treeHeight - this.depth;
    }
  }
  // Arrays to store Vertex objects for the two trees
  let vertices1 = [];
  let vertices2 = [];
  // Use a depth first search of each tree to fill Vertex arrays
  let dfs1 = tree1.depthFirstSearch({
     // Root to start the search from
    root:'#'+root1.id(),
    // Handler function called upon visiting each node: v is the current node, u the previous node, e their
    // connecting edge, i the index indicating the ith visited node and depth the depth of this node
    // We need the ids of the node, its parent as well as its depth to construct our Vertex object
    visit: function(v, e, u, i, depth) { 
      let parentId = (depth === 0) ? undefined : u.id();
      vertices1.push(new Vertex(v.id(), depth, parentId));
    }
  });
  let dfs2 = tree2.depthFirstSearch({
    root:'#'+root2.id(),
    visit: function(v, e, u, i, depth) {
      let parentId = (depth === 0) ? undefined : u.id();
      vertices2.push(new Vertex(v.id(), depth, parentId));
    }
  });
  // Calculate height of trees: This is simply the depth of the node of maximum depth (reduce returns an object)
  let tree1Height = vertices1.reduce( (a,b) => (b.depth > a.depth) ? b : a).depth;
  let tree2Height = vertices2.reduce( (a,b) => (b.depth > a.depth) ? b : a).depth;
  // Trees of different heights are certainly not isomorphic
  if (tree1Height !== tree2Height) {
    // console.log('Rooted trees have different heights and so not isomorphic.');
    return false;
  }
  // Javascript way of creating an array of empty arrays of fixed length treeHeight + 1
  let tree1Levels = Array.apply(null, Array(tree1Height + 1)).map( () => []);
  let tree2Levels = Array.apply(null, Array(tree1Height + 1)).map( () => []);
  // The ith entry of tree1Levels holds all the Vertex objects at level i in the first tree
  // Note we must set the vertex.treeHeight property as this is not done in the first depth searches
  vertices1.forEach( vertex => {vertex.treeHeight = tree1Height; tree1Levels[vertex.level].push(vertex);});
  vertices2.forEach( vertex => {vertex.treeHeight = tree1Height; tree2Levels[vertex.level].push(vertex);});
  // Initialisation complete, now to perform the algorithm, starting at level 0 and working upwards
  for (let i=0; i < tree1Levels.length - 1; i++) {
    // Assign integers according to .assignedTuple of each node (or 0 if leaf)
    assignLevelInts(tree1Levels[i]);
    assignLevelInts(tree2Levels[i]);
    // Trees must have the same number of vertices at each level in order to be isomorphic
    if (tree1Levels[i].length !== tree2Levels[i].length) {
      // console.log('Rooted trees do not have the same number of vertices at level ' + i + ' and so are not isomorphic.');
      return false;
    }
    // Flag to toggle if a reason trees are not isomorphic are found (why use flag? Cannot return from within a .forEach)
    let isomorphic = true;
    // Step 6. of algorithm - if two sequences of tuples are not identical, then halt; trees are not isomorphic.
    tree1Levels[i].forEach( (vertex, index) => {
      // compareVertexTuples(v1,v2) returns 0 if v1.assignedTuple and v2.assignedTuple are identical
      // We can compare Vertices at corresponding positions in tree1Levels & tree2Levels because these arrays have been sorted
      if (compareVertexTuples(vertex, tree2Levels[i][index]) !== 0) {
        // console.log('Sequences of tuples at level ' + i + ' do not match; rooted trees not isomorphic.');
        isomorphic = false;
      }
    });
    if (!isomorphic) {
      return false;
    }
    // Trees are isomorphic 'up to level i'; before moving to next level we need to populate the tuples of the branches
    // on the next level according to the value of assigned integers in the current level (Step 3. of algorithm)
    assignNextLevelTuples(tree1Levels[i], tree1Levels[i+1]);
    assignNextLevelTuples(tree2Levels[i], tree2Levels[i+1]);
  }
  // If ordered tuples are identical at level (treeHeight-1), then we are done; the trees are isomorphic (the roots of
  // each tree will be assigned the same integer in step 7.).
  return true;

  // Helper functions
  // Function to set .assignedInt for each Vertex in levelVerts array
  function assignLevelInts (levelVerts) {
    // Sort the tuple of each Vertex in non-decreasing order. See N.B. in assignNextLevelTuples below
    levelVerts.map( vertex => vertex.assignedTuple.sort() );
    // Now sort the VERTICES in levelVerts according to their TUPLES; non-decreasing variable length lexicographic sort
    levelVerts = levelVerts.sort(compareVertexTuples);
    // All leaves with .assignedTuple = [] (now at the start of levelVerts due to the sort) are assigned 0
    let intToAssign = 0;
    if (levelVerts[0].assignedTuple.length !== 0) {
      // There are no leaves as first Vertex in levelVerts has .assignedTuple !== []
      intToAssign++;
    }
    // Assigned integer to first Vertex manually (in below for we compare each Vertex's tuple to the PREVIOUS Vertex tuple)
    levelVerts[0].assignedInt = intToAssign;
    // Vertices with the same tupleare assigned the same integer
    // So assign value intToAssign until we find a difference from one tuple to the next; then intToAssign should be increment
    for (let j=1; j<levelVerts.length; j++) {
      // compareVertexTuples returns 0 only if tuples associated with Vertices are identical
      if (compareVertexTuples(levelVerts[j], levelVerts[j-1]) !== 0) {
        intToAssign++;
      }
      // Remember javascript passes objects by reference, so by editing levelVerts[j] we are changing the same Vertex
      // as it occurs in vertices1 or vertices2
      levelVerts[j].assignedInt = intToAssign;
    }
  }
  // Function to assign tuples to Vertices in nextLevelVerts using integers assigned to Vertices in levelVerts - step 3. of algorithm
  // This is done by iterating through levelVerts, and for each vertex a) finding its parent in nextLevelVerts and b) pushing
  // its .assignedInt onto .assignedTuple of this parent. N.B. This results in the tuples being unsorted, whilst the algorithm calls
  // for them to be sorted. We could sort the Vertices in levelVerts here according to their assigned integers (as suggested by Aho et. al.)
  // or, instead, just sort each tuple after assignment, as done above in assignLevelInts.
  function assignNextLevelTuples (levelVerts, nextLevelVerts) {
    levelVerts.forEach( vertex => {
      // find returns Vertex object in nextLEvelVerts with .id matching .parentId of current Vertex
      let father = nextLevelVerts.find( possibleParent => possibleParent.id === vertex.parentId );
      if (father === undefined) {
        console.error('Unable to find parent to vertex ' + vertex.id);
      }
      father.assignedTuple.push(vertex.assignedInt);
    });
  }
  // Function to compare two tuples of integers for 'variable length lexicographic sorting'
  // See Aho, Hopecroft, Jeffrey | Algorithm 3.2 p.80
  function compareVertexTuples (v1, v2) {
    // Tuples to be compared
    let tuple1 = v1.assignedTuple;
    let tuple2 = v2.assignedTuple;
    // Length of shortest tuple
    let minLength = (tuple1.length > tuple2.length) ? tuple2.length : tuple1.length;
    // Element=wise comparison of first minLength elements
    for (let j=0; j<minLength; j++) {
      if (tuple1[j] > tuple2[j]) {
        // tuple1 should occur AFTER tuple2 i.e. tuple2 takes precedence in a sort
        return 1;
      }
      if (tuple1[j] < tuple2[j]) {
        // tuple1 should occur BEFORE tuple2 i.e. tuple2 takes precedence in a sort
        return -1;
      }
    }
    // First minLength elements of tuples are identical, in which case the shortest tuple should occur first...
    if (tuple1.length > tuple2.length) {
      return 1;
    } else if (tuple1.length < tuple2.length) {
      return -1;
    }
    // Or tuples are the same length and hence identical (0 the sentinel value for equality of 'no sort required')
    return 0;
  }
}