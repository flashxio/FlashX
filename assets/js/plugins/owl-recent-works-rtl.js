var OwlRecentWorks = function () {

    return {

        ////Owl Recent Works v1
        initOwlRecentWorksV1: function () {
            jQuery(".owl-recent-works-v1").owlCarousel({
                loop: true,
                rtl: true,
                //nav: true,
                responsive: {
                    0:{
                        items: 1
                    },
                    600:{
                        items: 2
                    },
                    900:{
                        items: 3
                    },
                    1000:{
                        items: 4
                    }
                }
            });

            // Custom Navigation Events
            jQuery(".next-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('next.owl.carousel');
            })
            jQuery(".prev-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('prev.owl.carousel');
            })
        },  

        ////Owl Recent Works v2
        initOwlRecentWorksV2: function () {
            jQuery(".owl-recent-works-v1").owlCarousel({
                loop: true,
                rtl: true,
                //nav: true,
                responsive: {
                    0:{
                        items: 1
                    },
                    600:{
                        items: 2
                    },
                    1000:{
                        items: 3
                    }
                }
            });

            // Custom Navigation Events
            jQuery(".next-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('next.owl.carousel');
            })
            jQuery(".prev-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('prev.owl.carousel');
            })
        },        

        ////Owl Recent Works v3
        initOwlRecentWorksV3: function () {
            jQuery(".owl-recent-works-v1").owlCarousel({
                loop: true,
                rtl: true,
                //nav: true,
                responsive: {
                    0:{
                        items: 1
                    },
                    600:{
                        items: 2
                    },
                    1000:{
                        items: 1
                    }
                }
            });

            // Custom Navigation Events
            jQuery(".next-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('next.owl.carousel');
            })
            jQuery(".prev-v2").click(function(){
                jQuery(".owl-recent-works-v1").trigger('prev.owl.carousel');
            })
        },

    };
    
}();