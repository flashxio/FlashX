var CommentForm = function () {

    return {
        
        //Comment Form
        initCommentForm: function () {
    	    // Validation
	        $("#sky-form4").validate({
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
	                url:
	                {
	                    url: true
	                },
	                comment:
	                {
	                    required: true
	                },
	                captcha:
	                {
	                    required: true,
	                    remote: 'assets/plugins/sky-forms/version-2.0.1/captcha/process.php'
	                }
	            },
	                                
	            // Messages for form validation
	            messages:
	            {
	                name:
	                {
	                    required: 'Enter your name',
	                },
	                email:
	                {
	                    required: 'Enter your email address',
	                    email: 'Enter a VALID email'
	                },
	                url:
	                {
	                    email: 'Enter a VALID url'
	                },
	                comment:
	                {
	                    required: 'Please enter your comment'
	                },
	                captcha:
	                {
	                    required: 'Please enter characters',
	                    remote: 'Correct captcha is required'
	                }
	            },
	                                
	            // Ajax form submition                  
	            submitHandler: function(form)
	            {
	                $(form).ajaxSubmit(
	                {
	                    beforeSend: function()
	                    {
	                        $('#sky-form4 button[type="submit"]').attr('disabled', true);
	                    },
	                    success: function()
	                    {
	                        $("#sky-form4").addClass('submited');
	                    }
	                });
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