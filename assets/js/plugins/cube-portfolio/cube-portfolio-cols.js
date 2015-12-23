var CubePortfolio = function () {

    return {

        //Cube Portfolio 3 Columns
        initCubePortfolio3Col: function () {
        	jQuery(document).ready( function() {
			    jQuery('#grid-container').cubeportfolio({
			        // options
			        gridContainer.cubeportfolio({
				        mediaQueries: [{
				            width: 800,
				            cols: 3
				        }, {
				            width: 500,
				            cols: 2
				        }, {
				            width: 320,
				            cols: 1
				        }]
				    });
			    });
			});
        }
        
    };
    
}();

