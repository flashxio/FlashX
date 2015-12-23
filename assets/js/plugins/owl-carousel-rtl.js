var OwlCarousel = function () {

    return {

        //Owl Carousel
        initOwlCarousel: function () {
	        //Owl Slider v1
			var owl = jQuery(".owl-slider").owlCarousel({
                loop: true,
                rtl: true,
                margin: 30,
                responsive: {
                    0:{
                        items: 1
                    },
                    479:{
                        items: 2
                    },
                    600:{
                        items: 3
                    },
                    900:{
                        items: 4
                    },
                    1000:{
                        items: 5
                    }
                }
            });
            jQuery(".next-v1").click(function(){
                owl.trigger('next.owl.carousel');
            })
            jQuery(".prev-v1").click(function(){
                owl.trigger('prev.owl.carousel');
            })


	        //Owl Slider v2
			var owl2 = jQuery(".owl-slider-v2").owlCarousel({
                loop: true,
                rtl: true,
                margin: 30,
                responsive: {
                    0:{
                        items: 1
                    },
                    479:{
                        items: 2
                    },
                    600:{
                        items: 3
                    },
                    900:{
                        items: 4
                    },
                    1000:{
                        items: 5
                    }
                }
            });
            jQuery(".next-v2").click(function(){
                owl2.trigger('next.owl.carousel');
            })
            jQuery(".prev-v2").click(function(){
                owl2.trigger('prev.owl.carousel');
            })


	        //Owl Slider v3
			jQuery(".owl-slider-v3").owlCarousel({
                loop: true,
                rtl: true,
                responsive: {
                    0:{
                        items: 1
                    },
                    479:{
                        items: 2
                    },
                    600:{
                        items: 3
                    },
                    900:{
                        items: 4
                    },
                    1000:{
                        items: 5
                    },
                    1100:{
                        items: 9
                    }
                }
            });


	        //Owl Slider v4
			jQuery(".owl-slider-v4").owlCarousel({
                loop: true,
                rtl: true,
                nav: false,
                dots: true,
                dotsClass: "owl-pagination",
                dotClass: "owl-page",
                responsive: {
                    0:{
                        items: 1
                    },
                    479:{
                        items: 2
                    },
                    600:{
                        items: 3
                    }
                }
            });
			 

            //Owl Slider v5
            /*jQuery(document).ready(function() {
            var owl5 = jQuery(".owl-slider-v5");
                owl5.owlCarousel({
                    items:1,
                    itemsDesktop : [1000,1],
                    itemsDesktopSmall : [900,1],
                    itemsTablet: [600,1],
                    itemsMobile : [479,1]
                });
            });*/


            //Owl Slider v6
            jQuery(".owl-slider-v6").owlCarousel({
                loop: true,
                rtl: true,
                nav: false,
                dots: true,
                dotsClass: "owl-pagination",
                dotClass: "owl-page",
                responsive: {
                    0:{
                        items: 1
                    },
                    479:{
                        items: 2
                    },
                    979:{
                        items: 3
                    },
                    1000:{
                        items: 4
                    },
                    1100:{
                        items: 5
                    }
                }
            });

            //Owl Twitter v1
            /*jQuery(".owl-twitter-v1").owlCarousel({
                singleItem : true,
                slideSpeed : 1000,
                autoPlay : 10000,              
            });*/


            //Owl Testimonials v1
            jQuery(".owl-ts-v1").owlCarousel({
                loop: true,
                rtl: true,
                margin: 10,
                controlsClass: "owl-buttons",
                responsive: {
                    0:{
                        items: 1
                    }
                },
                navText: [,],
                nav: true,
            });


            //Owl Clients v1
            jQuery(".owl-clients-v1").owlCarousel({
                loop: true,
                rtl: true,
                merge: true,
                margin: 35,
                responsive: {
                    0:{
                        items: 1
                    },
                    300:{
                        items: 2
                    },
                    600:{
                        items: 3
                    },
                    900:{
                        items: 4
                    },
                    1000:{
                        items: 5
                    },
                    1280:{
                        items: 7
                    }
                }
            });


            //Owl Clients v2
            jQuery(".owl-clients-v2").owlCarousel({
                loop: true,
                rtl: true,
                margin: 10,
                responsive: {
                    0:{
                        items: 1
                    },
                    300:{
                        items: 2
                    },
                    600:{
                        items: 3
                    },
                    900:{
                        items: 4
                    },
                    1000:{
                        items: 5
                    }
                }
            });

            
            //Owl Video
            jQuery(".owl-video").owlCarousel({
                loop: true,
                rtl: true,
                dots: true,
                dotsClass: "owl-pagination",
                dotClass: "owl-page",
                responsive: {
                    0:{
                        items: 1
                    }
                }
            });            
		}
        
    };
    
}();