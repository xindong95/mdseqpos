/**
 * this script is responsible for displaying the motif table view
 */


function MotifTableView(motifs, container) {
    this.default_motifs = motifs; //this never changes
    this.motifs = motifs;
    this.container = container;
    this.colClickCntlr = new ColumnClickCntlr();
    this.thMap = {};
    this.imgMap = {};

    var outer = this;
    
    //to help the user restore the original view
    this.resetMtfs = function() {
	outer.setMotifs(outer.default_motifs);
    }

    this.setImg = function(field) {
	//Reset the other fields
	for (var i = 0; i < outer.fields; i++) {
	    outer.imgMap[outer.fields[i]].src = "down.png";
	}
	
	var imgSrc = 
	(outer.colClickCntlr.ascending) ? "down-red.png":"up-red.png";
	outer.imgMap[field].src = imgSrc;
    }

    this.fields = ["id", "factors", "consensus", "hits", "cutoff", 
		   "zscore", "pval", "position"];    
    this.makeHTML = function() {
	var tbl = document.createElement('table');
	tbl.className = "datatable"
	var header = document.createElement('tr');
	tbl.appendChild(header);

	for (var i = 0; i < outer.fields.length; i++) {
	    var field = outer.fields[i];
	    var tmp = document.createElement('th');
	    var span = document.createElement('span');
	    span.innerHTML = (field != "pval") ? field:"-10*log(pval)";
	    tmp.appendChild(span);
	    var img = document.createElement("img");
	    img.src = 'down.png';
	    tmp.appendChild(img);
	    tmp.onclick=function(f) { //have to curry it
		return function(event){
		    outer.colClickCntlr.setColumn(f);
		    outer.setImg(f);
		}
	    }(field);
	    header.appendChild(tmp);
	    outer.thMap[field] = tmp;
	    outer.imgMap[field] = img;
	}

	for (var i = 0; i < outer.motifs.length; i++) {
	    var newTr = document.createElement('tr');
	    newTr.onclick = function(n) { 
		return function(event) { 
		    motifModel.setCurrentMotif(outer.motifs[n]);}}(i);
	    //Alternating rows effect
	    newTr.className = (i % 2 == 0) ? "" : "altrow";
	    tbl.appendChild(newTr);
	    //associate the motif with the TR
	    //NOTE: must hash by string, not by object.
	    motifTRMap[outer.motifs[i].getid()] = newTr;
	    for (var j = 0; j < outer.fields.length; j++) {
		var tmp = document.createElement('td');
                if (outer.fields[j] == 'pval') { //SPECIAL CASE
                    //pvals are transformed
                    var pval = outer.motifs[i][outer.fields[j]];
                    tmp.innerHTML = -10*Math.log(pval);
                } else {
                    tmp.innerHTML = outer.motifs[i][outer.fields[j]];
		}
		newTr.appendChild(tmp);
	    }
	    //highlight the current selection -- BUT don't scroll b/c we 
	    //are column sorting
	    highlightMtfRow(false);
	}
	//put in the new calendar
	if (outer.container.childNodes[0] == null) {
	    outer.container.appendChild(tbl);
	} else {
	    outer.container.replaceChild(tbl, outer.container.childNodes[0]);
	}
	
    }
    
    this.setMotifs = function(newMotifs) {
	outer.motifs = newMotifs;
	outer.makeHTML();
    }

    this.init = function() {
	if (outer.motifs.length > 0) {
            outer.makeHTML();
	} else { //No motifs found, remove the reset btn
            var reset_btn = document.getElementById("reset_btn");
            reset_btn.parentNode.removeChild(reset_btn);
	}
    }
}

function ColumnClickCntlr() {
    this.lastClick = null;
    this.ascending = true;

    var outer = this;

    this.setColumn = function(column) {
	if (column == outer.lastClick) { //toggle
	    outer.ascending = !outer.ascending;
	} else {
	    outer.ascending = true;
	    outer.lastClick = column;
	}
	//make the call
	motifModel.sortMotifList(outer.lastClick, outer.ascending);
    }
}


//EDIT this --like the styling rules are NOT ok--well maybe
//NOTE: scroll parameter is now ignored b/c we dropped tree nav; OBSOLETE
function highlightMtfRow(scroll) {
    var prevMtf = motifModel.previousMotif;
    var currMtf = motifModel.currentMotif;

    if (prevMtf != null) {
	//unset it
	prevRow = motifTRMap[prevMtf.getid()];
	prevRow.style.border = "none";
    }
    
    if (currMtf != null) {
	//set the new row
	var currRow = motifTRMap[currMtf.getid()];
	currRow.style.border = "2px solid #CC6";
	/* OBSOLETE
	if (scroll) {
	    //automatically center the row by scrolling--window.scrollTo(x,y) 
	    //NOTE: we will only be changing the y coordinate--using the ratio
	    //of where the motif is found
	    var motifList = motifTableView.motifs;
	    var found = false;
	    var i = 0;
	    for (; !found && (i < motifList.length); i++) {
		if (motifList[i].id == currMtf.id) {
		    found = true;
		}
	    }
	    if (found) {
		//scroll using the RATIO of where it is in the list
		//window.scrollMaxY is the max scroll height of the document
		scrollToY = 
		    Math.floor((i / motifList.length) * window.scrollMaxY);
		window.scrollTo(0, scrollToY);
	    }
	}
	*/
    }
}

/******GLOBAL CALLS HERE ****/
var motifModel = parent.motifModel;

function motifListLstnr() {  
    motifTableView.setMotifs(motifModel.motifList);
}
motifModel.motifListEvent.register(motifListLstnr);

function currentMotifLstnr() {
    //highlight the current motif and scroll
    highlightMtfRow(true); 
}
motifModel.currentMotifEvent.register(currentMotifLstnr);

var motifTRMap = {};

var motifTableView = null;
var rst_btn = null; 
var species_menu = null;

function checkStrInList(s, list) {
    if (list != null) {
	for (var i = 0; i < list.length; i++) {
	    if (s == list[i]) {
		return true;
	    }
	}
    }
    return false;
}

function initPage(){
    //problem--this is undefined.
    motifTableView = new MotifTableView(motifModel.motifList, 
					document.getElementById('motif_table'));
    motifTableView.init();

    //A reset button to restore the table to a default state
    rst_btn = document.getElementById('reset_btn');
    if (rst_btn != null) {
	rst_btn.onclick = function(event) {
	    motifTableView.resetMtfs();
	}
    }

    //setup the species menu
    species_menu = document.getElementById('species_menu');
    if (species_menu != null) {
	//generate a list of unique species in the list
	species_list = ['Any'];
	if (motifModel.motifList != null) {
	    for (var i = 0; i < motifModel.motifList.length; i++) {
		species = motifModel.motifList[i].getspecies();
		for (var j = 0; j < species.length; j++) {
		    if (!checkStrInList(species[j], species_list)) {
			species_list.push(species[j]);
		    }
		}
	    }
	}
	//populate the menu with the species
	for (var i = 0; i < species_list.length; i++) {
	    new_opt = document.createElement("option");
	    new_opt.value = species_list[i];
	    new_opt.innerHTML = species_list[i];
	    species_menu.appendChild(new_opt);
	}
	
	//SET the onchange event
	species_menu.onchange = function(event) {
	    filterSpecies(this.options[this.selectedIndex].value);
	}
	
	//enable the menu
	species_menu.disabled = false;
    }
}

//Filters the motifModel according to the species specified
function filterSpecies(species) {
    //special case: 'Any' = no filtering
    if (species == "Any") {
	motifModel.setMotifList(motifModel.originalMotifList);
    } else {
	list = motifModel.originalMotifList;
	tmp = []
	for (var i = 0; i < list.length; i++) {
	    if (checkStrInList(species, list[i].species)) {
		tmp.push(list[i]);
	    }
	}
	motifModel.setMotifList(tmp);
    }
}
