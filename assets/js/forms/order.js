var OrderForm = function () {

    return {
        
        //Order Form
        initOrderForm: function () {
	        // Datepickers
	        $('#start').datepicker({
	            dateFormat: 'dd.mm.yy',
	            prevText: '<i class="icon-chevron-left"></i>',
	            nextText: '<i class="icon-chevron-right"></i>',
	            onSelect: function( selectedDate )
	            {
	                $('#finish').datepicker('option', 'minDate', selectedDate);
	            }
	        });
	        $('#finish').datepicker({
	            dateFormat: 'dd.mm.yy',
	            prevText: '<i class="icon-chevron-left"></i>',
	            nextText: '<i class="icon-chevron-right"></i>',
	            onSelect: function( selectedDate )
	            {
	                $('#start').datepicker('option', 'maxDate', selectedDate);
	            }
	        });
	        
	        // Validation
	        $("#sky-form1").validate({
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
	                phone:
	                {
	                    required: true
	                },
	                interested:
	                {
	                    required: true
	                },
	                budget:
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
	                phone:
	                {
	                    required: 'Please enter your phone number'
	                },
	                interested:
	                {
	                    required: 'Please select interested service'
	                },
	                budget:
	                {
	                    required: 'Please select your budget'
	                }
	            },

	            // Ajax form submition
	            submitHandler: function(form)
	            {
	                $(form).ajaxSubmit(
	                {
	                    beforeSend: function()
	                    {
	                        $('#sky-form1 button[type="submit"]').addClass('button-uploading').attr('disabled', true);
	                    },
	                uploadProgress: function(event, position, total, percentComplete)
	                {
	                    $("#sky-form1 .progress").text(percentComplete + '%');
	                },
	                    success: function()
	                    {
	                        $("#sky-form1").addClass('submited');
	                        $('#sky-form1 button[type="submit"]').removeClass('button-uploading').attr('disabled', false);
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