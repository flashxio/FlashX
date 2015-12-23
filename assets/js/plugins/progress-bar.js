// Created by: Farukham: (https://github.com/farukham/Bootstrap-Animated-Progress-Bars)

var ProgressBar = function () {

    return {

    	//Progress Bar Horizontal
    	initProgressBarHorizontal: function () {
	        jQuery(document).ready(function () {
			    jQuery('.progress').each(function () {
			        jQuery(this).appear(function() {
			          	jQuery(this).animate({opacity:1,left:"0px"},800);
			          	var b = jQuery(this).find(".progress-bar").attr("data-width");
			          	jQuery(this).find(".progress-bar").animate({
			            	width: b + "%"
			          	}, 100, "linear");
			        }); 
			    });   
			});
	    },

	    //Progress Bar Vertical
    	initProgressBarVertical: function () {
	        jQuery(document).ready(function () {
			    jQuery('.progress').each(function () {
			        jQuery(this).appear(function() {
			          	jQuery(this).animate({opacity:1,left:"0px"},800);
			          	var b = jQuery(this).find(".progress-bar").attr("data-height");
			          	jQuery(this).find(".progress-bar").animate({
			            	height: b + "%"
			          	}, 100, "linear");
			        }); 
			    });   
			});
	    }
	
	};
	
}();    