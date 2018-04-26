// November 18, 2014
// For use with plot_clone_size_distribution_figures.py

// Change these ALL-CAPS variables to customize the plot --------------

var x_max = 30.0; // default to false; set to false to use all the data
var y_max = 10000000; // 1e7;

var suppress_obs_zero = false; // if actual data, don't
					   // show the circle that
					   // coresponds to n0
					   // (unobserved clones)

// Some general/miscellaneous variables -------------------------------

var y_min = 7e-1;

var dimension_base = 500; // overall size scale of the plot

var filename = "/Users/ramy/Downloads/Recon-master/test_sample_1_plotfile.txt";

var dataset; // into which we'll put the data; declare it globally so
	     // if you do anything outside of d3.tsv, data is available

var neg = function(d) { if (d < 1) { return "-" } else { return "" } };

var epsilon = 1e-10;

// Define sizes of elements in terms of dimension_base ---------------

var marker_size = dimension_base/50;

var stroke_width = {
    "axis": dimension_base/150,
    "marker": dimension_base/200
};

var margin = { // see http://bl.ocks.org/mbostock/3019563
    top: dimension_base/10, 
    right: marker_size + 2*stroke_width.marker,
    bottom: dimension_base/3.75,
    left: dimension_base/3
};

var width = dimension_base - margin.right - margin.left;
var height = 0.9*dimension_base - margin.top - margin.bottom;

var numXticks = Math.round(2 + width/200);
var numYticks = Math.round(2 + height/300);
// numMinorXticks is defined in d3.tsv because we need the data loaded

var font_size = {
    "ticks": dimension_base/17,
    "axis": dimension_base/13,
    "crosses": dimension_base/13 // in px
};

var exponent_ratio = 0.5;
var exponent_size = exponent_ratio*font_size.axis;
var exponent_displacement = 1.1*(exponent_ratio - 1);

var font_weight = function(d) {
    var divisor = 2.2;
    if (dimension_base/divisor < 100) { return 100 }
    else if (dimension_base/divisor > 900) { return 900 }
    else { return dimension_base/divisor }
}


// Append the svg and load the data ----------------------------------

var svg = d3.select("body")
    .append("svg")
    .attr({
	    "width": width + margin.left + margin.right,
	    "height": height + margin.top + margin.bottom
	})
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    ;

d3.tsv(filename, function(error, data) {
	dataset = data;

	// filter out y-values below the min; saves on empty x-axis
	dataset = data.filter(function(d) {
		if (d.fit >= y_min) { return d; }
	    });

	// find the max x value, as a default
	var x_max_data = d3.max(dataset, function(d) { return +d.clone_size; } );
	if (Boolean(x_max) == false) { x_max = x_max_data; };

	// filter out data with x-values above the max
	dataset = data.filter(function(d) {
		if (d.clone_size <= x_max) { return d; }
	    });

	// find the max y value, as a default
	dataset = dataset.filter(function(d) {
		if (d.clone_size <= 48) { return d; }
	    });

	// if no y_max is given, find the max y value, as a default
	if (y_max == false) { 
	    var y_max_data = d3.max(dataset, function(d) { return +d.fit; } );
	    if (suppress_obs_zero == "true") {
		y_max = dataset[0].upper_limit;
	    }
	    else {
		y_max = d3.max(dataset, function(d) { return +d.fit; } );
	    };
	    y_max = Math.log(y_max) / Math.LN10;
	    y_max = Math.ceil(y_max);
	    y_max = Math.pow(10, y_max);	
	};

	var numMinorXticks;
	if (x_max <= 30) { numMinorXticks = x_max; }
	else { numMinorXticks = 2*numXticks; };

	// Scale data to plot size ----------------------------------
	
	var axis_gap_divisor = 25 // the smaller the divisor, the bigger the gap between the origins of the axes

	var x_scale = d3.scale.linear()
	    .domain ( [0, x_max] )
	    .range([dimension_base/axis_gap_divisor, width])
	    // .nice()
	    ;

	var y_scale = d3.scale.log()
	    // .domain([1e-4, d3.max(dataset, function(d) { return +d.fit; } ) ])
	    // fit will almost always be bigger because of n0; should max(d.fit, d.sample) for safety, though
	    .domain([1e0, y_max])
	    .range([height-dimension_base/axis_gap_divisor, 0])
	    ;
	
	// Define and plot axes -------------------------------------

	// x-axis
	var xAxis = d3.svg.axis()
	    .scale(x_scale)
	    .orient("bottom")
	    .ticks(numXticks)
	    .tickSize(dimension_base/70)
	    .tickPadding(dimension_base/20)
	    ;
	
	svg.append("g")
	    .attr({
		    "class": "axis",
		    "transform": "translate(0," + height + ")"
		    })
	    .style({
		    "stroke-width": stroke_width.axis,
		    "font-size": font_size.ticks + "px"
		    })
	    .call(xAxis)
	    .append("text")
	    .attr({
		    "x": width/2,
		    "y": margin.bottom - font_size.axis/4,
		    "text-anchor": "middle"
		  })
	    .style("font-size", font_size.axis + "px")
	    .text("clone size, sample")
	    ;

	// a second x-axis is used to provide minor tickmarks
	var xAxisMinorTicks = d3.svg.axis()
	    .scale(x_scale)
	    .orient("bottom")
	    .ticks(numMinorXticks)
	    .tickSize(dimension_base/70)
	    .tickPadding(dimension_base/20)
	    ;

	svg.append("g")
	    .attr({
		    "class": "axis",
		    "transform": "translate(0," + height + ")"
		    })
	    .style({
		    "stroke-width": stroke_width.axis,
		    "font-size": font_size.ticks + "px",
		    "fill": "none"
		    })
	    .call(xAxisMinorTicks)
	    ;

	// y-axis
	var yAxisTicks = [];
	for ( var i = 0; i < (Math.log(y_max) / Math.LN10) + 1; i++ ) {
	    yAxisTicks.push(Math.pow(10, i));
	};
	
	var yAxis = d3.svg.axis()
	    .scale(y_scale)
	    .orient("left")
	    .tickValues(yAxisTicks)//.slice(0, numYticks))
	    .tickSize(dimension_base/70)
	    .tickPadding(dimension_base/30)
	    ;
	
	svg.append("g")
	    .attr("class", "axis")
	    .style({
		    "stroke-width": stroke_width.axis,
		    "font-size": font_size.ticks + "px"
		    })
	    .call(yAxis)
	    .selectAll(".tick text")
	    .text(null)
	    .filter(powerOfTen)
	    .text(10)
	    .append("tspan")
	    .attr("dy", exponent_displacement + "em")
	    .style("font-size", exponent_size + "px")
	    .text(function(d) { return Math.round(Math.log(d) / Math.LN10); })
	    ;

	function powerOfTen(d) {
	    return d / Math.pow(10, Math.ceil(Math.log(d) / Math.LN10 - 1e-12)) === 1;
	}

	// Error bars (included only if suppress_obs_zero is on) -------------

	console.log(+dataset[0].upper_limit);
	console.log(isNaN(+dataset[0].upper_limit));
	console.log(typeof(+dataset[0].upper_limit));
	if (!isNaN(+dataset[0].upper_limit) ) {

	    var error_bars = svg.append("g") // error_bars is a group
		.attr("class", "error_bars")
		;

	    error_bars.append("line") // Lower limit
		.attr({
			"x1": x_scale(0)-marker_size,
			"y1": y_scale(+dataset[0].lower_limit),
			"x2": x_scale(0)+marker_size,
			"y2": y_scale(+dataset[0].lower_limit),
			"stroke": "black",
			"stroke-width": stroke_width.marker
			})
		;
	    error_bars.append("line") // Upper limit
		.attr({
			"x1": x_scale(0)-marker_size,
			"y1": y_scale(+dataset[0].upper_limit),
			"x2": x_scale(0)+marker_size,
			"y2": y_scale(+dataset[0].upper_limit),
			"stroke": "black",
			"stroke-width": stroke_width.marker
			})
		;
	}
	;
	    
	// Plot data ------------------------------------------------------

	var circles = svg.selectAll("circles") // ...which don't exist yet
	    .data(dataset) // for each thing in dataset
	    .enter()
	    .append("circle")
	    ;

	circles.attr({

		cx: function(d) { 
		    if (suppress_obs_zero == true) {
			if (+d.clone_size != 0) { return x_scale(+d.clone_size); }
			else { return x_scale(-100000); } // poor man's way of supporessing data; should do it when reading in data
		    }
		    else { return x_scale(+d.clone_size); }
		},

		// the following if... else catches zeroes
		cy: function(d) { 
		    if (+d.sample != 0) { return y_scale(+d.sample); }
		    else if (y_scale(+d.sample) < y_min) { return y_scale(epsilon); }
		    else { return y_scale(epsilon); }
		},
		r: marker_size
		})
	    .style({
		    "stroke-width": stroke_width.marker,
		    "fill": "none",
		    "stroke": function(d) {
			if (d.clone_size == 0) { return "#f00"; } // make n0 red
			else { return "#000"; }
		    }
		})
	    ;

	var crosses = svg.selectAll("text.cross") // selects text of
						  // class "cross"
	    .data(dataset)
	    .enter()
	    .append("text")
	    .append("tspan")
	    .attr("class", "cross")
	    .style({
		    "font-size": font_size.crosses + "px",
		    "font-weight": font_weight
		})
	    ;

	crosses.text("+")
	    .attr({
		    rotate: 45,
			x: function(d) { return x_scale(+d.clone_size) - marker_size/1.5 + stroke_width.marker/2; },
		    y: function(d) { 
			if (+d.fit >= y_min) { return y_scale(+d.fit); } // used to include stroke_width.marker/2
			else { return y_scale(epsilon); }
		    },
		    "text-anchor": "middle"
		})
	    ;

	// add the error bars for n0
	

    });

