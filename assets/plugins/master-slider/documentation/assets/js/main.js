$(document).ready(function(){
	
	SyntaxHighlighter.all();
	
	CTween.fadeOut($('#loading') , 400 , true);
	CTween.fadeIn($('#main-container').css('display' , 'block') , 400);
	
	$(".easedemo").each(function(index){
				var $this = $(this);
				var tween;
				var pbtn = $('<div class="playbtn"></div>').appendTo($this).click(function(){
					if(tween) tween.reset();
					box.css('left' , 0);
					tween = CTween.animate(box , 1800 , {left:440} , {ease:$this.data('ease')});
				});
				var box = $('<div class="easebox"></div>').appendTo($this);
			});
	
	// init pages
	
	var pages = $('.content-section');
	var content = $('#content').addClass('float');
	pages.css('display','none').addClass('float');
	
	var currentPage;
	
	function updatePage(){
		var hash = window.location.hash;
		if(hash === '' || hash === undefined) hash = '#intro';
		var page = $(hash);
		showPage(page);
		if(currentPage)	hidePage(currentPage);
		currentPage = page;
		content.scrollTop(0);
	}
	
	var hide_tween ;
	
	var showPage =  function(page){
		$('a[href=#'+page.attr('id')+']').parent().addClass('active');
		
		page.css({opacity:0 , display:'' , position:'relative' , top:'0px' , left:'0px'});
		CTween.setPos(page,{x:500});
		CTween.animate(page , 500 , {opacity:'1' , x:0} , {ease:'easeOutQuart'});
		
		if(hide_tween && ((hide_tween.$element && hide_tween.$element[0] === page[0]) || (hide_tween[0] === page[0]))){
			hide_tween.stop(true);
		}
		
	};
	
	var hidePage = function(page){
		$('a[href=#'+page.attr('id')+']').parent().removeClass('active');
		page.css({position:'absolute'});
		var to = {opacity:'0' , x : 500};
		hide_tween = CTween.animate(page , 500 , to , {ease:'easeOutQuart', complete:function(){
			page.css('display' , 'none');
		}});
	};
	

	$(window).on('hashchange' , updatePage);
	updatePage();
	
	$(window).on('resize' , onresize);
	
	function onresize(){
		content.height($(window).height() - $('.header').height() - parseInt($('#content').css('padding-top'))*2 - 1);
	}
	onresize();
	
	
	var sidebar = $('#sidebar'),
		side_w	= $('.toc').outerWidth(),
		show	= true,
		togg	= $('.toggle'), 
		base_marg = parseInt(content.css('margin-left'));
		
	togg.on('click',function(){
		togg.toggleClass('out');
		if(show){
			show = false;
			CTween.animate(sidebar , 500 , {left:-side_w+'px'} , {ease:'easeOutQuart'});
			CTween.animate(content , 500 , {marginLeft:base_marg-side_w} , {ease:'easeOutQuart'});
		}else{
			show = true;
			CTween.animate(sidebar , 500 , {left:'0px'} , {ease:'easeOutQuart'});
			CTween.animate(content , 500 , {marginLeft:base_marg} , {ease:'easeOutQuart'});
		}
		
	});
	
});