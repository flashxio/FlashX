<!doctype html>
<html>
<head>
   <meta charset="UTF-8">
   <link rel="shortcut icon" href="../.favicon.ico">
   <title>Cliff-Linux Contents</title>
   <style>
     *
     	{padding:0;
	 margin:0;
	 -o-transition:all 200ms linear;
	 -moz-transition:all 200ms linear;
	 -webkit-transition:all 200ms linear;
	 transition:all 200ms linear;}

     html,body 
     	{color:#333;
	 font-family: "Lucida Console", Courier, monospace;
	 font-size:14px;}
 
     #container
	{margin:0 auto;
	 width:900px;
	 margin-top:20px;
	 padding-top:10px;
	 border:1px solid #EEE;
	 border-radius:10px;
	 position:relative;}
     
     table 
	{background-color:#F3F3F3;
	 border-collapse:collapse;
	 width:100%;
	 margin:15px 0;}

     th
     	{text-align:left;
	 background-color:#FE4902;
	 color:#FFF;
	 font-weight:bold;
	 padding:7px 0 5px 22px;
	 font-size:14px;
	 letter-spacing:1px;
	 font-family: Verdana, Arial, Helvetica, sans-serif;
	 margin-bottom:5px;}
     
     th span
         {font-size:9px; 
	  letter-spacing:0;}
     td
	{padding:0px 0;
	 position:relative;
	 overflow:hidden;}

     a
     	{text-decoration:none;}
	
     table a
     	{padding: 6px 0 6px 24px;
	 color:#663300;
	 display:block;
	 width:100%;
	 height:100%;}

     table tr td:first-of-type a
     	{background: url(./.images/file.png) no-repeat left 50%;}

     table tr td:not(:first-of-type) a
     	{background-image:none !important;} 

     tr:nth-of-type(odd) 
         {background-color:#E6E6E6;}

     tr:hover td
         {background-color:#CACACA;}

     tr:hover a 
         {color:#000;}

     h1
	{font-size:18px;
	 font-weight:bold;
	 padding:0 0 10px 10px;
	 text-align:center;}
     
     h2
     	{font-size:16px;
     	 text-align:center;}

     /* icons for file types (icons by famfamfam.) */

     /* images */
     table tr td:first-of-type a[href$=".jpg"], 
     table tr td:first-of-type a[href$=".png"], 
     table tr td:first-of-type a[href$=".gif"], 
     table tr td:first-of-type a[href$=".svg"], 
     table tr td:first-of-type a[href$=".jpeg"]
      		{background-image: url(./.images/image.png);}

     /* pdfs */
     table tr td:first-of-type a[href$=".pdf"] 
     		{background-image: url(./.images/pdf.gif);}
     
     /* zips */
     table tr td:first-of-type a[href$=".zip"] 
     		{background-image: url(./.images/zip.png);}

     /* css */
     table tr td:first-of-type a[href$=".css"] 
     		{background-image: url(./.images/css.png);}

     /* docs */
     table tr td:first-of-type a[href$=".doc"],
     table tr td:first-of-type a[href$=".docx"],
     table tr td:first-of-type a[href$=".ppt"],
     table tr td:first-of-type a[href$=".pptx"],
     table tr td:first-of-type a[href$=".pps"],
     table tr td:first-of-type a[href$=".ppsx"],
     table tr td:first-of-type a[href$=".xls"],
     table tr td:first-of-type a[href$=".xlsx"]
     		{background-image: url(./.images/office.png)}
      
     /* videos */
     table tr td:first-of-type a[href$=".avi"], 
     table tr td:first-of-type a[href$=".wmv"], 
     table tr td:first-of-type a[href$=".mp4"], 
     table tr td:first-of-type a[href$=".mov"], 
     table tr td:first-of-type a[href$=".m4a"]
     		{background-image: url(./.images/video.png);}

     /* audio */
     table tr td:first-of-type a[href$=".mp3"], 
     table tr td:first-of-type a[href$=".ogg"], 
     table tr td:first-of-type a[href$=".aac"], 
     table tr td:first-of-type a[href$=".wma"] 
     		{background-image: url(./.images/audio.png);}

     /* web pages */
     table tr td:first-of-type a[href$=".html"],
     table tr td:first-of-type a[href$=".xml"]
     		{background-image: url(./.images/xml.png);}
      
     table tr td:first-of-type a[href$=".php"] 
     		{background-image: url(./.images/php.png);}

     table tr td:first-of-type a[href$=".js"] 
     		{background-image: url(./.images/script.png);}

     /* directories */
     table tr.dir td:first-of-type a
     		{background-image: url(./.images/folder.png);}

   </style>
</head>

<body>
<div id="container">
	<h1>Cliff-Linux Contents (192.168.1.74)</h1> 
	<h2>Your IP Address is: <?php echo getenv('REMOTE_ADDR'); ?></h2>
	
   	<table>
		<tr>
			<th>Filename</th>
			<th>Type</th>
			<th>Size <span>(bytes)</span></th>
			<th>Date Modified</th>
		</tr>
	<?php
	 // Opens directory
	 $myDirectory=opendir(".");

	 // Gets each entry
	 while($entryName=readdir($myDirectory)) {
	   $dirArray[]=$entryName;
	 }

	 // Finds extensions of files
	 function findexts ($filename)
	 {
	   $filename=strtolower($filename);
	   $exts=split("[/\\.]", $filename);
	   $n=count($exts)-1;
	   $exts=$exts[$n];
	   return $exts;
	 }

	 // Closes directory
	 closedir($myDirectory);

	 // Counts elements in array
	 $indexCount=count($dirArray);

	 // Sorts files
	 sort($dirArray);

	 // Loops through the array of files
	 for($index=0; $index < $indexCount; $index++) {

	 // Allows ./?hidden to show hidden files
		if($_SERVER['QUERY_STRING']=="hidden")
		{$hide="";
		 $ahref="./";
		 $atext="Hide";}
		else
		{$hide=".";
		 $ahref="./?hidden";
		 $atext="Show";}
	       if(substr("$dirArray[$index]", 0, 1) != $hide) {

	 // Gets File Names
	       $name=$dirArray[$index];
	       $namehref=$dirArray[$index];

	 // Gets Extensions 
	       $extn=findexts($dirArray[$index]); 

 	 // Gets file size 
	       $size=number_format(filesize($dirArray[$index]));

	 // Gets Date Modified Data
	       $modtime=date("M j Y g:i A", filemtime($dirArray[$index]));

	 // Prettifies File Types
		switch ($extn){
			case "png": $extn="PNG Image"; break;
			case "jpg": $extn="JPG Image"; break;
			case "svg": $extn="SVG Image"; break;
			case "gif": $extn="GIF Image"; break;
			case "ico": $extn="Windows Icon"; break;

			case "txt": $extn="Text File"; break;
			case "htm": $extn="HTML File"; break;
			case "php": $extn="PHP Script"; break;
			case "js": $extn="Javascript"; break;
			case "css": $extn="Stylesheet"; break;
			case "pdf": $extn="PDF Document"; break;

			case "zip": $extn="ZIP Archive"; break;

			default: $extn=strtoupper($extn)." File"; break;
		}

	// Separates directories
	       if(is_dir($dirArray[$index]))
	       {
		       $extn="&lt;directory&gt;"; 
		       $size="&lt;directory&gt;"; 
		       $class="dir";
		}
		else 
		{
			$class="file";
		}

	// Cleans up . and .. directories 
		if($name=="."){$name=". (Current Directory)";}
		if($name==".."){$name=".. (Parent Directory)";}

	 // Print 'em
	 print("
	       	<tr class='$class'>
			<td><a href='$namehref'>$name</a></td>
			<td><a href='$namehref'>$extn</a></td>
			<td><a href='$namehref'>$size</a></td>
			<td><a href='$namehref'>$modtime</a></td>
		</tr>");
	   }
	 }
	?>

	</table>

	<h2><?php print("<a href='$ahref'>$atext hidden files</a>"); ?></h2>
</div>
</body>
</html>