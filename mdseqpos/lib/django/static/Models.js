/*
 * Models.js
 *
 * Copyright (c) 2007 Len Taing
 *
 * Last Modified: Time-stamp: <2009-05-01 17:49:44 Tao Liu>
 *
 * Description:
 * Model classes hold all of the data that a page requires.  The data is 
 * accessed through "getter" functions, e.g. Foo.getName().  The data is 
 * updated through "setter" functions, e.g. Foo.setName("bob").  When the data
 * changes, the model notifies its listeners of the update.  These listeners
 * are usually "Views" that may refresh themselves based on the changes in the
 * data model.
 *
 * ModelEvent is the class that handles all of the model event notifications.
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
    
    this.unregister = function(listenerFn) {
	var i = outer.listeners.indexOf(listenerFn);
	if (i != -1) { outer.listeners.splice(i, 1);}
    }
    
    this.notify = function() {
	outer.listeners.each(function(lstnr) {lstnr(outer.modelRef);});
    }
}

/**
 * Class: ModelFactory
 * Description: This Class Factory generates Model classes.  Model classes 
 * have three responsibilities: 1. store data, 2. have set/get fns to access 
 * the data, and 3. have ModelEvents that are invoked when the data changes
 *
 * This factory takes two arrays as arguments: getsetters are fields that we 
 * we want define getter AND setter functions for (they are the "read/write"
 * fields).  They also have associated ModelEvents--when the setter function
 * changes the field's value, the event is invoked.
 *
 * getters are "readonly" i.e. they only have associated "getter" fns.
 *
 * NOTE: field names MUST follow java variable name conventions: firstword
 * lowercase, every other word after is upper-cased.
 *
 * generates methods like: (assuming the field is called "firstField") 
 * getFirstField, setFirstField, firstFieldEvent.
 *
 * The factory returns a CLASS that expects one argument: an object that 
 * associates each field name w/ a the given value, e.g. {foo:5, bar:true},
 * OR EntryModel(entryId, placeId,...) --> EntryModel({entryId:,placeId:});
 *
 * @param getsetters - instance vars we want to be read&writable
 * @param getters - instance vars we want to be readonly
 * @return a new ModelClass
 *
 */
function ModelFactory(getsetters, getters) {
    return function(paramObj) {
	var outer = this;
	this.fields = getsetters.concat(getters);
	
	//create the getter functions
	for (var i = 0; i < getters.length; i++) {
	    this[getters[i]] = paramObj[getters[i]];
	    var name = upperCase1stLtr(getters[i]);
	    //we need a closure to grab the REF to the field value
	    this["get"+name] = function(field) { 
		return function(){ return outer[field]; }
	    }(getters[i]);

	}

	//set the incoming params: e.g. this.entryId = paramObj.entryId
	//and create get functions
	for (var i = 0; i < getsetters.length; i++) {
	    this[getsetters[i]] = paramObj[getsetters[i]];
	    var name = upperCase1stLtr(getsetters[i]);
	    //Again we need a closure for help
	    this["get"+name] = function(field) {
		return function(){return outer[field]; }
	    }(getsetters[i]);
	    
	    //create the event and set functions
	    this[getsetters[i]+"Event"] = new ModelEvent(this);
	    this["set"+name] = function(field) {
		return function(param) {
		    if (param != outer[field]) {
			outer[field] = param;
			outer[field+"Event"].notify();
		    }
		}
	    }(getsetters[i])
	}
    }
}

/**
 * Fn: takes a string and returns the string w/ the first letter uppercased
 * NOTE: maybe I should move the following to a file called StringUtils.
 */
function upperCase1stLtr(s) {
    if (s == null) { return null;}
    if (s.length > 0) {
	var ltr = s.substring(0,1);
	var rest = s.substring(1,s.length);
	return ltr.toUpperCase() + rest;
    } else {//empty string
	return "";
    }
}
