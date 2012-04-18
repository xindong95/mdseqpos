/*
 * focus frame: this script is responsible for creating the view in the 
 * focus frame.
 */

/******GLOBAL CALLS HERE ****/
//does this work??? YES: no need for an onload.
//PROBLEM: the onload call mucks up this timing--the frame is loaded before
//the var is init.
//SOLN: drop the init.
var motifModel = parent.motifModel;

function currMtfLstnr() {
    var container = document.getElementById('focus');
    displayMotif(motifModel.currentMotif, container, true, true, document);
}

motifModel.currentMotifEvent.register(currMtfLstnr);

function init() {
    //motifModel.currentMotifEvent.notify();
    if (motifModel.motifList.length == 0) { //no motifs found
        var container = document.getElementById('focus');
        displayNoMotifs(container, document);
    }
}

//Notifies the user that no motifs were found
function displayNoMotifs(container, dom) {
    var newH2 = dom.createElement('h2');
    newH2.innerHTML = "Sorry no motifs were found.  <br/>Please try re-running the tool <br/> with wider search parameters";
    if (container.childNodes[0] == null) {
        container.appendChild(newH2);
    } else {
        container.replaceChild(newH2, container.childNodes[0]);
    }
}

function displayMotif(motif, container, showNewWindowBtn, showPssmBtn, dom) {
    var newTable = dom.createElement('table');
    var newTr = dom.createElement('tr');
    newTable.appendChild(newTr);

    //Add the LOGO
    var logoTd = dom.createElement("td");
    newTr.appendChild(logoTd);
    var img = dom.createElement("img");
    logoTd.appendChild(img);
    img.src = motif.getlogoImg();    
    img.className="logoimg";

    //display the info
    var newTd = dom.createElement('td');
    newTr.appendChild(newTd);
    var newUl = dom.createElement('ul');
    newTd.appendChild(newUl);

    var fields = ["id", "factors", "consensus", "hits", "cutoff",
                  "zscore", "pval", "position"];
    for (var i = 0; i < fields.length; i++) {
        var newLi = dom.createElement('li');
	var field = fields[i];
        var fieldname = (field != "pval") ? field : "pval";
	newLi.innerHTML = fieldname+":"+motif[field];
        newUl.appendChild(newLi);
    }

    if (showNewWindowBtn) {
	var btn = dom.createElement('input');
	btn.type = "button";
	btn.value = "show motif in new window";
	btn.onclick = function(event) { 
	    popUpWindow(motif);
	}
	newUl.appendChild(btn);
    }

    if (showPssmBtn) {
	//add a 'show pssm btn'
	var pssmBtn = dom.createElement('input');
	pssmBtn.type = "button";
	pssmBtn.value = "show pssm in a new window";
	pssmBtn.onclick = function(event) {
	    popUpPssm(motif);
	}
	newUl.appendChild(pssmBtn);
    }

    if (container.childNodes[0] == null) {
	container.appendChild(newTable);
    } else {
	container.replaceChild(newTable, container.childNodes[0]);
    }
}

function popUpPssm(motif) {
    var win = window.open("focus_frame_2.html", '',"width=300, height=600,"+
			  "resizable=yes, scrollbars=yes, toolbar=no,"+
			  "location=no, menubar=no, status=yes");
    //NEED to create the PSSM matrix AFTER the page is loaded
    win.onload = function(event) {
	win.document.title = motif.id + " PSSM";
	//load the pssm_popup.css style
	var css = win.document.createElement('link');
	css.href = 'pssm_popup.css';
	css.rel = "stylesheet";
	css.type = "text/css";
	win.document.body.appendChild(css);

	//BUILD the table
	var tbl = win.document.createElement('table');
	tbl.className = "datatable";
	win.document.body.appendChild(tbl);

	var header = win.document.createElement('tr');
	tbl.appendChild(header);
	//BUILD the header
	var order = ['', 'A', 'C', 'G', 'T'];
	for (var i = 0; i < order.length; i++) {
	    var th = win.document.createElement('th');
	    th.innerHTML = order[i];
	    header.appendChild(th);
	}
	var pssm = "["; //represent the pssm as a matrix
	for (var i = 0; i < motif.pssm.length; i++) {
	    var tr = win.document.createElement('tr');
	    tbl.appendChild(tr);
	    //The nucleotide index
	    var tmp = win.document.createElement('td');
	    tmp.innerHTML = i+1;
	    if (i % 2 == 0) { tr.className="altrow"; }
	    tr.appendChild(tmp);
	    var row_pssm = "["; //represent the pssm row
	    for (var j = 0; j < motif.pssm[i].length; j++) {
		var td = win.document.createElement('td');
		tr.appendChild(td);
		td.innerHTML = motif.pssm[i][j];
		
		//build row_pssm
		row_pssm += (j != 0) ? ", " : "";
		row_pssm += motif.pssm[i][j];
	    }
	    row_pssm += "]"; //close row
	    pssm += (i != 0) ? ",\n" : "";
	    pssm += row_pssm;
	}
	pssm += "]"; //close the matrix
	
	//var build the txt dump
	var div = win.document.createElement('div');
	win.document.body.appendChild(div);
	var p = win.document.createElement('p');
	p.innerHTML="If you would like to input this motif into our \'Screen all motifs\' tool: <br/>OPTION A: <br/>1. copy the text below, and <br/>2. paste it into \"PSSM Raw Text\" in the Motif:Screen All Motifs tool<br/>--OR--<br/></br>OPTION B:<br/>1. copy and paste the txt below into a text file <br/>2. upload the file onto cistrome, and <br/>3. run Motif:Screen All Motifs selecting the text file as the PSSM file";
	p.style.fontSize="10pt";
	div.appendChild(p);
	var pre = win.document.createElement('pre');
	pre.innerHTML = pssm;
	pre.style.fontSize="10pt";
	div.appendChild(pre);	
    }
}

function popUpWindow(motif) {
    var win = window.open("focus_frame_2.html", '',"width=600, height=180,"+
			  "resizable=yes, scrollbars=no,toolbar=no,"+
			  "location=no, menubar=no, status=yes");

    //NEED to create the motif elements AFTER the page load
    win.onload = function(event) {
	//BUILD the new window's DOM
	win.document.title = motif.id;
	var tmp = win.document.createElement('div');
	win.document.body.appendChild(tmp);
	displayMotif(motif, tmp, false, false, win.document);
    }
}
