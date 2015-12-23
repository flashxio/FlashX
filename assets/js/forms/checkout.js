var CheckoutForm = function () {

    return {
        
        //Checkout Form
        initCheckoutForm: function () {
	        // Masking
	        $('#cvv').mask('999', {placeholder:'X'});
	        $('#card').mask('9999-9999-9999-9999', {placeholder:'X'});
	        $('#year').mask('2099', {placeholder:'X'});
	        
	        // Add validation method
	        $.validator.addMethod("creditcard", function(value, element, param) {
	            if( /[^0-9\-]+/.test(value) )
	            {
	                return false;
	            }
	        },
	        $.validator.format('Please enter a valid credit card number.'));
	        
	        // Validation
	        $('#sky-form').validate({
	            // Rules for form validation
	            rules:
	            {
	                fname:
	                {
	                    required: true
	                },
	                lname:
	                {
	                    required: true
	                },
	                email:
	                {
	                    required: true,
	                    email: true
	                },
	                phone:
	                {
	                    required: true
	                },
	                country:
	                {
	                    required: true
	                },
	                city:
	                {
	                    required: true
	                },
	                code:
	                {
	                    required: true,
	                    digits: true
	                },
	                address:
	                {
	                    required: true
	                },
	                name:
	                {
	                    required: true
	                },
	                card:
	                {
	                    required: true,
	                    creditcard: true
	                },
	                cvv:
	                {
	                    required: true,
	                    digits: true
	                },
	                month:
	                {
	                    required: true
	                },
	                year:
	                {
	                    required: true,
	                    digits: true
	                }
	            },
	                                
	            // Messages for form validation
	            messages:
	            {
	                fname:
	                {
	                    required: 'Please enter your first name'
	                },
	                lname:
	                {
	                    required: 'Please enter your last name'
	                },
	                email:
	                {
	                    required: 'Please enter your email address',
	                    email: 'Please enter a VALID email address'
	                },
	                phone:
	                {
	                    required: 'Please enter your phone number'
	                },
	                country:
	                {
	                    required: 'Please select your country'
	                },
	                city:
	                {
	                    required: 'Please enter your city'
	                },
	                code:
	                {
	                    required: 'Please enter code',
	                    digits: 'Digits only please'
	                },
	                address:
	                {
	                    required: 'Please enter your full address'
	                },
	                name:
	                {
	                    required: 'Please enter name on your card'
	                },
	                card:
	                {
	                    required: 'Please enter your card number'
	                },
	                cvv:
	                {
	                    required: 'Enter CVV2',
	                    digits: 'Digits only'
	                },
	                month:
	                {
	                    required: 'Select month'
	                },
	                year:
	                {
	                    required: 'Enter year',
	                    digits: 'Digits only please'
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