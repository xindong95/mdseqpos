/******GLOBAL CALLS HERE ****/
var motifModel = parent.motifModel;

var st;
function init() {
    var json = motifModel.getMotifTree();
    //Create a new canvas instance.
    var canvas = new Canvas('mycanvas', {
	    //Where to inject canvas. Any HTML container will do.
	    'injectInto':'infovis',
	    //Set width and height, default's to 200.
	    'width': 600,
	    'height': 180,
	    //Set a background color in case the browser
	    //does not support clearing a specific area.
	    'backgroundColor': '#222'
	});
    //Create a new ST instance
    st= new ST(canvas, {
	    //Set node and edge colors
	    //Set overridable=true to be able
	    //to override these properties
	    //individually
	    constrained: false,

	    Node: {
		overridable: true,
		color: '#ccb'
	    },
	    Edge: {
		overridable: true,
		color: '#ccb'
	    },
	    //Add an event handler to the node when creating it.
	    onCreateLabel: function(label, node) {
		label.id = node.id;
		label.innerHTML = node.name;

		//FOR the LEAF NODES, the label should be the id
		//NOTE: the 'children' field is replaced by 'adjacencies'
		var ct = 0;
		for (prop in node.adjacencies) { ct++; }
		if (ct == 1) { //ITS a leaf node b/c only adjacent to its parnt
		    label.innerHTML = node.id;
		}
		    
		label.onclick = function(){
		    st.onClick(node.id);
		    //ONCLICK, set it as the current motif.
		    motifModel.setCurrentMotif(node.data.motif);
		};
	    },
	    //This method is called right before plotting
	    //a node. It's useful for changing an individual node
	    //style properties before plotting it.
	    //The data properties prefixed with a dollar
	    //sign will override the global node style properties.
	    onBeforePlotNode: function(node) {
		//add some color to the nodes in the path between the
		//root node and the selected node.
		if (node.selected) {
		    node.data.$color = "#ff7";
		} else {
		    delete node.data.$color;
		}
	    },

	    //This method is called right before plotting
	    //an edge. It's useful for changing an individual edge
	    //style properties before plotting it.
	    //Edge data properties prefixed with a dollar sign will
	    //override the Edge global style properties.
	    onBeforePlotLine: function(adj){
		if (adj.nodeFrom.selected && adj.nodeTo.selected) {
		    adj.data.$color = "#eed";
		    adj.data.$lineWidth = 3;
		}
		else {
		    delete adj.data.$color;
		    delete adj.data.$lineWidth;
		}
	    }
	});
    //load json data
    st.loadJSON(json);
    //compute node positions and layout
    st.compute();
    //optional: make a translation of the tree
    //Tree.Geometry.translate(st.tree,
    //			    new Complex(-200, 0), "startPos");
    //Emulate a click on the root node.
    st.onClick(st.root);

    //Make the root image clickable
    makeRootClickable();
}

/* listens to the change in the current node, and redisplays the tree
 */
function highlightTreeNode() {
    //NOTE: the st id is set to the motif id
    var currMotif = motifModel.currentMotif;
    //WARNING:global call
    st.onClick(currMotif.getid());
}

/* each time the user clicks, this function with draw a representation of 
the root to node path; it listens for current Motif events
*/
function drawNavBar() {
    var currMotif = motifModel.currentMotif;
    var pathToRoot = [];
    var currNode = currMotif.treeNode;
    //To build the path to root, walk up the parent links until parent=null
    while (currNode.parent != null) {
	//alert(currNode.parent);
	pathToRoot.push(currNode);
	currNode = currNode.parent.treeNode;
    }
    var rootToNode = pathToRoot.reverse();
    
    var navbar = document.getElementById('nav_bar');
    //clear the span
    navbar.innerHTML = "";
    //build these up with images
    for (var i = 0; i < rootToNode.length; i++) {
	//append the connector image
	var tmp = document.createElement('img');
	tmp.src = "cnnctor.png";
	navbar.appendChild(tmp);

	//then the node image
	tmp = document.createElement('img');
	tmp.src = "node.png";
	tmp.onclick = function(node) {
	    return function(event) {
		motifModel.setCurrentMotif(node.node);
	    }
	}(rootToNode[i])
	navbar.appendChild(tmp);
    }
}

//This function makes the root image clickable
function makeRootClickable() {
    var rootImg = document.getElementById('root_img');
    rootImg.onclick = function(node) {
	return function(event) {
	    motifModel.setCurrentMotif(node.node);
	}
    }(motifModel.motifTree);
}


motifModel.currentMotifEvent.register(highlightTreeNode);
motifModel.currentMotifEvent.register(drawNavBar);

