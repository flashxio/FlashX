

$(document).ready(function() {

	// CSS transitions
	$('body').addClass('anim');

	// Init syntax highlighter
	SyntaxHighlighter.defaults['toolbar'] = false;
	SyntaxHighlighter.all();

	// jump to deeplink
	$('a[data-goto]').click(function(){
		$('#sidebar li[data-deeplink="'+$(this).data('goto')+'"]').click();
	});

	// Init sidebar

	$('#sidebar li').click(function(e){

		var subNav = '';
		var hash = '';

		$(this).addClass('active').siblings().removeClass('active');

		if( $(this).parent().closest('li').length ){

			e.stopPropagation();

			var pli = $(this).parent().closest('li');

			// checking if the parent main menuitem is active or not

			if( !pli.hasClass('active') ){

				pli.addClass('active').siblings().removeClass('active');
				$('#content > section').removeClass('active').eq( pli.index() ).addClass('active');
			}

			// submenu

			$('#content > section.active > article').removeClass('active').eq( $(this).index() ).addClass('active');

			// setting variables

			hash = $(this).data('deeplink');
			subNav = ' / ' + $('#sidebar .active .active').text();

		}else{

			// main menu

			$('#content > section').removeClass('active').eq( $(this).index() ).addClass('active');

			var dl = $(this);

			// selecting the first submenu

			if( $(this).find('li').length ){

				$(this).find('li').removeClass('active');
				dl = $(this).find('li:first').addClass('active');
				$('#content > section.active > article').removeClass('active').eq(0).addClass('active');

				// setting variable

				subNav = ' / ' + $('#sidebar .active .active').text();
			}

			// setting variable

			hash = dl.data('deeplink');
		}

		// changing curnav

		$('#curnav').text( $('#sidebar li.active > span').text() + subNav );

		// changing hash

		document.location.hash = hash;

		// fading in

		$('#content > section.active').css('display','none').fadeIn(200);
	});

	if( document.location.hash ){

		$('#sidebar li[data-deeplink="'+document.location.hash.split('#')[1]+'"]').click();
	}

	if( typeof layerSliderTransitions !== 'undefined' ){

		// Init transition gallery
		transitionGallery.init();
	}
});


var transitionGallery = {

	init : function() {

		var self =  this;

		// Add transition list
		self.appendTransitions(layerSliderTransitions['t2d'], $('#slide-transitions-2d tbody'));
		self.appendTransitions(layerSliderTransitions['t3d'], $('#slide-transitions-3d tbody'));

		// Show transitions
		jQuery('.ls-transition-list').on('mouseenter', 'a', function() {
			self.showTransition(this);
		});

		// Hide transitions
		jQuery('.ls-transition-list').on('mouseleave', 'a', function() {
			self.hideTransition(this);
		});
	},

	appendTransitions : function(transitions, target) {
		for(c = 0; c < transitions.length; c+=2) {

			// Append new table row
			var tr = jQuery('<tr>').appendTo(target).append('<td class="c light"></td><td></td><td class="c"></td><td></td>');


			// Append transition col 1 & 2
			tr.children().eq(0).text((c+1));
			tr.children().eq(1).append( jQuery('<a>', { 'href' : '#', 'html' : transitions[c]['name'], 'rel' : 'tr'+(c+1) } ) )
			if(transitions.length > (c+1)) {
				tr.children().eq(2).text((c+2));
				tr.children().eq(3).append( jQuery('<a>', { 'href' : '#', 'html' : transitions[(c+1)]['name'], 'rel' : 'tr'+(c+2) } ) )
			}
		}
	},


	showTransition : function(el) {

		// Get transition index
		var index = jQuery(el).attr('rel').substr(2)-1;

		// Create popup
		jQuery('body').prepend( jQuery('<div>', { 'class' : 'ls-popup' })
			.append( jQuery('<div>', { 'class' : 'inner ls-transition-preview' }))
		);

		// Get popup
		var popup = jQuery('.ls-popup');

		// Get viewport dimensions
		var v_w = jQuery(window).width();

		// Get element dimensions
		var e_w = jQuery(el).width();

		// Get element position
		var e_l = jQuery(el).offset().left;
		var e_t = jQuery(el).offset().top;

		// Get toolip dimensions
		var t_w = popup.outerWidth();
		var t_h = popup.outerHeight();

		// Position tooltip
		popup.css({ top : e_t - t_h - 60, left : e_l - (t_w - e_w) / 2  });

		// Fix top
		if(popup.offset().top < 20) {
			popup.css('top', e_t + 75);
		}

		// Fix left
		if(popup.offset().left < 20) {
			popup.css('left', 20);
		}

		// Get transition class
		var trclass = jQuery(el).closest('.ls-transition-list').attr('id');

		// Built-in 3D
		if(trclass == 'slide-transitions-3d') {
			var trtype = '3d';
			var trObj = layerSliderTransitions['t'+trtype+''][index];

		// Built-in 2D
		} else if(trclass == 'slide-transitions-2d') {
			var trtype = '2d';
			var trObj = layerSliderTransitions['t'+trtype+''][index];
		}

		// Init transition
		popup.find('.inner').lsTransitionPreview({
			transitionType : trtype,
			transitionObject : trObj,
			imgPath : 'assets/img/',
			skinsPath: '../layerslider/skins/',
			delay : 100
		});
	},

	hideTransition : function(el) {

		// Stop transition
		jQuery('.ls-popup').find('.inner').lsTransitionPreview('stop');

		// Remove transition
		jQuery('.ls-popup').remove();
	},
};
