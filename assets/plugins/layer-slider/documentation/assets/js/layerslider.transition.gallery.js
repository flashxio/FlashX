(function( $ ) {
	$.fn.lsTransitionPreview = function(options) {

		return this.each( function(){
			if( typeof options === 'string' ){

				var tpData = $(this).data('transitionPreview');

				switch( options ){
					case 'stop':
						tpData.stop();
					break;
				}
			}else{
				new transitionPreview(this, options);
			}
		});
	};

	var transitionPreview = function(el, options){

		var tp = this;
		tp.$el = $(el);
		tp.$el.data('transitionPreview', tp);

		tp.init = function(){

			// Parse settings

			var settings = $.extend({
				width: 300,
				height: 150,
				delay: 100,
				imgPath: '../assets/img/',
				skinPath: '../layerslider/skins/',
				transitionType: '2d',
				transitionObject: null
			}, options );

			// Add slider HTML markup

			$(el).append( $('<div>', { 'class' : 'transitionpreview', 'style' : 'width: '+settings.width+'px; height: '+settings.height+'px;'})
				.append( $('<div>', { 'class' : 'ls-layer', 'data-ls' : 'slidedelay: '+settings.delay+';'})
					.append( $('<img>', { 'src' : ''+settings.imgPath+'sample_slide_1.png', 'class' : 'ls-bg'})))
				.append( $('<div>', { 'class' : 'ls-layer', 'data-ls' : 'slidedelay: '+settings.delay+';'})
					.append( $('<img>', { 'src' : ''+settings.imgPath+'sample_slide_2.png', 'class' : 'ls-bg'})))
			);

			// Initialize the slider

			$(el).find('.transitionpreview').layerSlider({
				showCircleTimer : false,
				pauseOnHover : false,
				skin : 'noskin',
				slidedelay : 100,
				skinsPath : settings.skinsPath,
				slideTransition : {
					type : settings.transitionType,
					obj : settings.transitionObject
				}
			});
		};

		tp.stop = function(){

			$(el).find('.transitionpreview').layerSlider('forceStop').remove();
		};

		tp.init();
	};
}( jQuery ));
