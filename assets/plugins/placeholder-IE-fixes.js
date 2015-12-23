var PlaceholderIEFixes = function () {

    return {
        
        //Placeholder IE Fixes
        initPlaceholderIEFixes: function () {
	        if (jQuery.browser.msie && jQuery.browser.version.substr(0, 1) < 9) { // ie7&ie8

	            jQuery('input[placeholder], textarea[placeholder]').each(function () {
	                var input = jQuery(this);
	                var inputCloneTypePass = $('<input type="text">');
	                var displayCss = input.css('display');

	                if ( input.val() == '' ) {

	                        if ( input.attr('type') == 'password' ) {
	                        $.each(input.get(0).attributes, function(v,n) {
	                                n = n.nodeName||n.name;
	                                if ( n != 'type' && n != 'name' ) {
	                                v = input.attr(n); // relay on $.fn.attr, it makes some filtering and checks
	                                if(v != undefined && v !== false) {
	                                        inputCloneTypePass.attr(n,v);
	                                }
	                                }
	                        });

	                        input.css('display', 'none');
	                        inputCloneTypePass
	                                .appendTo(input.parent())
	                                .val(input.attr('placeholder'))
	                                .focus(function () {
	                                            if (inputCloneTypePass.val() == inputCloneTypePass.attr('placeholder')) {
	                                                inputCloneTypePass.css('display', 'none');
	                                                input.css('display', displayCss);
	                                                input.focus();
	                                            }
	                                        });
	                        }

	                        input.val(input.attr('placeholder'));
	                }

	                jQuery(input).focus(function () {
	                    if (input.val() == input.attr('placeholder')) {
	                        input.val('');
	                    }
	                });

	                jQuery(input).blur(function () {
	                    if (input.val() == '' || input.val() == input.attr('placeholder')) {
	                        if (input.attr('type') == 'password') {
	                                inputCloneTypePass.css('display', displayCss);
	                                input.css('display', 'none');
	                        }else {
	                                input.val(input.attr('placeholder'));
	                        }
	                    }
	                });
	            });
	        }
        }

    };

}();        