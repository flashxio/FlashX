var OwlRecentWorks = function () {

    return {

        ////Owl Recent Works v1
        initOwlRecentWorksV1: function () {
            jQuery(document).ready(function() {
            var owl = jQuery(".owl-recent-works-v1");
                owl.owlCarousel({
                    items: [4],
                    itemsDesktop : [1000,4],
                    itemsDesktopSmall : [900,3],
                    itemsTablet: [600,2],
                    itemsMobile : [479,1],
                    slideSpeed: 1000
                });

                // Custom Navigation Events
                jQuery(".next-v2").click(function(){
                    owl.trigger('owl.next');
                })
                jQuery(".prev-v2").click(function(){
                    owl.trigger('owl.prev');
                })
            });            
        },  

        ////Owl Recent Works v2
        initOwlRecentWorksV2: function () {
            jQuery(document).ready(function() {
            var owl = jQuery(".owl-recent-works-v1");
                owl.owlCarousel({
                    items: [3],
                    itemsDesktop : [1000,3],
                    itemsDesktopSmall : [900,2],
                    itemsTablet: [600,2],
                    itemsMobile : [479,1],
                    slideSpeed: 1000
                });

                // Custom Navigation Events
                jQuery(".next-v2").click(function(){
                    owl.trigger('owl.next');
                })
                jQuery(".prev-v2").click(function(){
                    owl.trigger('owl.prev');
                })
            });            
        },        

        ////Owl Recent Works v3
        initOwlRecentWorksV3: function () {
            jQuery(document).ready(function() {
            var owl = jQuery(".owl-recent-works-v1");
                owl.owlCarousel({
                    items: [1],
                    itemsDesktop : [1000,1],
                    itemsDesktopSmall : [900,2],
                    itemsTablet: [600,2],
                    itemsMobile : [479,1],
                    slideSpeed: 1000
                });

                // Custom Navigation Events
                jQuery(".next-v2").click(function(){
                    owl.trigger('owl.next');
                })
                jQuery(".prev-v2").click(function(){
                    owl.trigger('owl.prev');
                })
            });            
        }

    };
    
}();