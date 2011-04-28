/*
 * This file contains the motif model.  It holds and processes all of the 
 * data for the page.  The different views (of the model, i.e. the different 
 * frames of the output), will register with the model and listen for 
 * model changes
 *
 */

/**
 * Class: ModelEvent
 * Description: This class does two important things: 1. register listeners
 * for a model event, and 2. when the model event is "invoked" (via the 
 * notify method), all of the registered listeners are informed.
 *
 * @param: modelRef -- the object that is passed as the message of notify
 * NOTE: it is customary for modelRef to be MODEL that uses the Model Event
 */
function ModelEvent(modelRef) {
    this.modelRef = modelRef;
    this.listeners = [];
    var outer = this;
    
    this.register = function(listenerFn) {
	var i = outer.listeners.indexOf(listenerFn);
	//make sure there are no duplicates
	if (i == -1) { outer.listeners.push(listenerFn);}
    }
    
    //WARNING: i think this is dependent on protoype and shouldn't be used!
    this.unregister = function(listenerFn) {
	var i = outer.listeners.indexOf(listenerFn);
	if (i != -1) { outer.listeners.splice(i, 1);}
    }
    
    this.notify = function() {
	for (var i = 0; i < outer.listeners.length; i++) {
	    outer.listeners[i](outer.modelRef);
	}
	//outer.listeners.each(function(lstnr) {lstnr(outer.modelRef);});
    }
}

/**
 * Class: Motif
 * this function contains all of the information for any individual motif
 * fields: id, factors, consensus sequence, pssm - motif matrix, logoImg,
 * hit score, cutoff score, zscore, pvalue, and gene group
 * 
 * note: we don't have gene group information yet
 * 
 * The input should be a hashtable of the fields and their associated value
 */
//DELETE THIS???
function MotifHelper(id, factors, consensus, pssm, logoImg, hits, 
		     cutoff, zscore, pval) {
    return new Motif({"id":id, "factors":factors, "consensus":consensus,
		"pssm":pssm, "logoImg":logoImg, "hits":hits, "cutoff":cutoff,
		"zscore":zscore, "pval":pval});
}

//this function returns the maximum score of the pssm row
function max_score(row) {
    if ((row == null) || (row.length != 4)) {
	return 0;
    }
    
    var n1 = (row[0] > row[1]) ? row[0] : row[1];
    var n2 = (row[2] > row[3]) ? row[2] : row[3];
    return (n1 > n2) ? n1 : n2;
}

//A dictionary of nucleotides to IUPAC symbols
var IUPAC_DICT = {'A':'A', 'C':'C', 'G':'G', 'T':'T', //Singles     
                  'AC':'M', 'AG':'R', 'AT':'W', 'CG':'S', 'CT':'Y', 'GT':'K',
                  'ACG':'V', 'ACT':'H', 'AGT':'D', 'CGT':'B',
                  'ACGT':'N'};

//a function that translates a row in the pssm to its corresponding PSSM symbol
function IUPAC(row) {
    if (row == null) {
        return "*";
    } else if (row.length != 4) {
        return "*";
    }
    var key="";
    var order=["A","C","G","T"];
    for (var i = 0; i < row.length; i++) {
        if (row[i] > 0.24) { //IF it exceeds the cutoff, add it to the key
            key += order[i];
        }
    }
    return IUPAC_DICT[key];
}

//this function returns the shortened version of the consensus sequence
//it takes the top 5 nucleotides based on scores and represents the gaps
//in the sequence w/ the number of nucleotides omitted
function shortenConsensus(consensus, row_scores) {
    //construct an ordered list of scores and indexes
    var list = [];
    for (var i = 0; i < row_scores.length; i++) {
	var tmp = {'score':row_scores[i], 'index':i};
	//insert it in list- splice(index, number, item)
	//insertion sort
	var found = false;
	var j = 0;
	for (; j < list.length && !found; j++) {
	    found = tmp.score > list[j].score;
	}
	if (found) { //insert
	    list.splice(j - 1, 0, tmp);
	} else { //append
	    list.push(tmp);
	}
    }

    //now we take the top 5 indices
    var indices = [list[0].index, list[1].index, list[2].index, 
		   list[3].index, list[4].index];
    //SORT them -- bubble sort!
    for (var i = 0; i < indices.length; i++) {
	for (var j = 0; j < (indices.length - 1); j++) {
	    if (indices[j] > indices[j+1]) { //swap!
		var tmp = indices[j+1];
		indices[j+1] = indices[j];
		indices[j] = tmp;
	    }
	}
    }

    //construct the string:
    var last = -1;
    var s = "";
    var diff;
    for (var i = 0; i < indices.length; i++) {
	diff = indices[i] - last;
	if (diff > 1) {//print the number
	    s += diff - 1;
	}
	s += consensus.charAt(indices[i]);
	last = indices[i];
    }
    //add info about the end IF last is not the last nucleotide
    diff = (consensus.length - 1) - last;
    if (diff > 0) { s += diff;}

    return s;
}

function Motif(paramObjs) {	
    this.fields = ["id", "factors", "entrezs", "refseqs", "species", 
		   "consensus", "pssm", "logoImg", "hits", "cutoff", "zscore",
		   "pval", "position"];
    var outer = this;

    /* OBSOLETE
    this.treeNode = null; //this corresponds to the motif's treeNode object
    this.setTreeNode = function(node) {
    	outer.treeNode = node;
    }
    */

    //1. set the values
    //e.g. this is short-hand for: this.id = paramObjs.id
    //2. create "get" functions, e.g. getid()
    for (var i = 0; i < outer.fields.length; i++) {
	this[outer.fields[i]] = paramObjs[outer.fields[i]];
	this["get"+outer.fields[i]] = function(field) {
	    return function(){return outer[field]; }
	}(outer.fields[i]);
    }
    
    //construct the consensus string:                                      
    this.consensus = "";
    var row_scores = [];
    for (var r = 0; r < this.pssm.length; r++) {
        this.consensus += IUPAC(this.pssm[r]);
	row_scores.push(max_score(this.pssm[r]));
    }

    //IF the string exceeds the max-chars(12), then we must shorten
    if (this.consensus.length > 12) {
	this.consensus = shortenConsensus(this.consensus, row_scores);
    }

    this.toString = function() {
	tmp = "";
	for (var i = 0; i < outer.fields.length; i++) {
	    tmp += outer.fields[i]+":"+paramObjs[outer.fields[i]]+"\n";
	}
	return tmp;
    }
}

/**
 * Class: MotifModel
 * Description: this model contains the following information-
 * 1. an ordered list of motifs
 * 2. the current motif
 * 3. the motif tree - DROPPED b/c the tree nav window is gone
 *
 * note: the current motif is set to null, so no need to pass it in as 
 * a param to the constructor
 */

function MotifModel(motifList) {
    this.originalMotifList = motifList; //READ-ONLY!
    this.motifList = motifList;
    //this.motifTree = motifTree; OBSOLETE

    this.currentMotif = null;
    this.previousMotif = null;
    
    var outer = this;

    //events
    this.motifListEvent = new ModelEvent(this);
    //this.motifTreeEvent = new ModelEvent(this);
    this.currentMotifEvent = new ModelEvent(this);

    this.setMotifList = function(newMotifList) {
	if (outer.motifList != newMotifList) {
	    outer.motifList = newMotifList;
	    outer.motifListEvent.notify();
	}
    }

    this.setCurrentMotif = function(newMotif) {
	if (outer.currentMotif != newMotif) {
	    outer.previousMotif = outer.currentMotif;
	    outer.currentMotif = newMotif;
	}
	//WE always notify--even if it stays the same, so that the VIEWS
	//like, highlightMtfRow can update
	outer.currentMotifEvent.notify();
    }

    //sorts the motif list according to the field parameter, and 
    //sends out a notification to the listeners
    //ascending: true- ascending order, descending otherwise
    this.sortMotifList = function(field, ascending) {
	//remove any rows with the value None before sorting
	var nones = [];
	var tmp = [];
	for (var i = 0; i < outer.motifList.length; i++) {
            if ((outer.motifList[i][field] == 'None') ||
                (field =='factors' && outer.motifList[i][field].length == 0)) {
		nones.push(outer.motifList[i]);
	    } else {
		tmp.push(outer.motifList[i]);
	    }
	}
	
	tmp = quicksort(tmp, field);
	tmp = (ascending) ? tmp : tmp.reverse();
	outer.motifList = tmp.concat(nones);
	outer.motifListEvent.notify();
    }
    
    this.toString = function() {
	tmp = "";
	for (var i = 0; i < outer.motifList.length; i++) {
	    tmp += motifList[i].toString();
	}
	return tmp;
    }

    //returns the json tree that jit expects, i.e.
    // {id: .., name:.., data:.., children: [...]}
    /* OBSOLETE
    this.getMotifTree = function() { 
	return outer.treeBuilder(outer.motifTree); 
    }
    
    //given a node in the motifTree, it returns the json 
    //description of that subtree.
    this.treeBuilder = function(node) {
	if (node == null) {
	    return {};
	}
	var children = [];
	for (var i = 0; i < node.children.length; i++) {
	    children.push(outer.treeBuilder(node.children[i]));
	}
	return {"id":node.node.getid(),
		"name":node.node.getconsensus(),
		"data":{"motif":node.node, "parent":node.parent},
		"children": children};
    }
    */
}

function quicksort(motifList, field) {
    var less = [];
    var greater = [];
    if (motifList.length <= 1) {
	return motifList;
    }
    var pivot = Math.floor(Math.random()*motifList.length);
    for (var i = 0; i < motifList.length; i++) {
	//remove the pivot
	if (i != pivot) { //otherwise you can get infinite recur. w/ equal
	    //should i use the getter functions here?
	    //alert(motifList[i][field]+","+motifList[pivot][field]);
	    if (motifList[i][field] <= motifList[pivot][field]) {
		less.push(motifList[i]);
	    } else {
		greater.push(motifList[i]);
	    }
	}
    }
    var tmp1 = quicksort(less, field);
    var tmp2 = [motifList[pivot]];
    var tmp3 = quicksort(greater, field);
    return tmp1.concat(tmp2.concat(tmp3));
}

