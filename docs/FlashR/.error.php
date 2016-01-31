<?php

$status=$_SERVER['REDIRECT_STATUS'];
$codes=array(
       400 => array('400 Bad Request', 'The request cannot be fulfilled due to bad syntax.'),
       401 => array('401 Login Error', 'It appears that the password and/or user-name you entered was incorrect. <a href="#" onclick="window.location.reload()">Click here</a> to return to the login page.'),
       403 => array('403 Forbidden', 'The server has refused to fulfill your request.'),
       404 => array('404 Not Found', 'Whoops, sorry, but the document you requested was not found on this server.'),
       405 => array('405 Method Not Allowed', 'The method specified in the Request-Line is not allowed for the specified resource.'),
       408 => array('408 Request Timeout', 'Your browser failed to send a request in the time allowed by the server.'),
       414 => array('414 URL To Long', 'The URL you entered is longer than the maximum length.'),
       500 => array('500 Internal Server Error', 'The request was unsuccessful due to an unexpected condition encountered by the server.'),
       502 => array('502 Bad Gateway', 'The server received an invalid response from the upstream server while trying to fulfill the request.'),
       504 => array('504 Gateway Timeout', 'The upstream server failed to send a request in the time allowed by the server.'),
);
 
$errortitle = $codes[$status][0];
$message = $codes[$status][1];

?>

<!doctype html>
<html>
<head>
	<title>That's an Error!</title>
	<style>
	  html 
	{color:#333;
	 font-family: "Lucida Console", Courier, monospace;
	 font-size:14px;
	 background:#eeeeee;}
 
	.content
	{margin:0 auto;
	 width:1000px;
	 margin-top:20px;
	 padding:10px 0 10px 0;
	 border:1px solid #EEE;
     background: none repeat scroll 0 0 white;
     box-shadow: 0 5px 10px -5px rgba(0, 0, 0, 0.5);
     position: relative;
}

	h1
		{font-size:18px;
		 text-align:center;}

	h1.title 
		{color:red;}
	
	h2
		{font-size:16px;
		 text-align:center;}
	
	p 
		{text-align:center;}

	hr
		{border:#fe4902 solid 1px;}

	</style>
</head>

<body>

	<div class="content">
	<h1>Sorry, but that's an error!</h1>
	<h1 class="title"><?php echo $errortitle; ?></h1>
	<hr>
	<p><?php echo $message;?></p>
	</div>

</body>
</html>