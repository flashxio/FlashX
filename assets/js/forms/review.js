var ReviewForm = function () {

    return {
        
        //Review Form
        initReviewForm: function () {
	        // Validation
	        $("#sky-form2").validate({
	            // Rules for form validation
	            rules:
	            {
	                name:
	                {
	                    required: true
	                },
	                email:
	                {
	                    required: true,
	                    email: true
	                },
	                review:
	                {
	                    required: true,
	                    minlength: 20
	                },
	                quality:
	                {
	                    required: true
	                },
	                reliability:
	                {
	                    required: true
	                },
	                overall:
	                {
	                    required: true
	                }
	            },
	                                
	            // Messages for form validation
	            messages:
	            {
	                name:
	                {
	                    required: 'Please enter your name'
	                },
	                email:
	                {
	                    required: 'Please enter your email address',
	                    email: 'Please enter a VALID email address'
	                },
	                review:
	                {
	                    required: 'Please enter your review'
	                },
	                quality:
	                {
	                    required: 'Please rate quality of the product'
	                },
	                reliability:
	                {
	                    required: 'Please rate reliability of the product'
	                },
	                overall:
	                {
	                    required: 'Please rate the product'
	                }
	            },                  
	            
	            // Do not change code below
	            errorPlacement: function(error, element)
	            {
	                error.insertAfter(element.parent());
	            }
	        });
        }

    };

}();