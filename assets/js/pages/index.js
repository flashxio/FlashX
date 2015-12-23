var Index = function () {

    return {
        
        //Parallax Slider
        initParallaxSlider: function () {
			$('#da-slider').cslider({
			    current     : 0,    
			    // index of current slide
			     
			    bgincrement : 50,  
			    // increment the background position 
			    // (parallax effect) when sliding
			     
			    autoplay    : false,
			    // slideshow on / off
			     
			    interval    : 4000  
			    // time between transitions
			});
        },

        //Layer Slider
        initLayerSlider: function () {
		    $(document).ready(function(){
		        jQuery("#layerslider").layerSlider({
			        skin: 'fullwidth',
			        responsive : true,
			        responsiveUnder : 960,
			        layersContainer : 960,
			        skinsPath: 'assets/plugins/layer-slider/layerslider/skins/'
			    });
		    });     
        },

        //Revolution Slider
        initRevolutionSlider: function () {
		    var revapi;
	        jQuery(document).ready(function() {
	            revapi = jQuery('.tp-banner').revolution(
	            {
	                delay:9000,
	                startwidth:1170,
	                startheight:500,
	                hideThumbs:10
	            });
	        });
        }

    };
}();        