// Toggle a division as expanded or collapsed.
// Also toggle the arrow icon.
// Refer to the division and image by their IDs.
//
// "Collapsed" material is hidden using the
// display property in CSS.

// Used by adaptproduct function (see below)
// to support adaptive doc in the Windows
// version of the Help Browser.
var adaptiveIds = new Array();
function toggleexpander(blockid, arrowid) {
    arrow = document.getElementById(arrowid);
    block = document.getElementById(blockid);
    if (block.style.display == "none") {
        // Currently collapsed, so expand it.
        block.style.display = "block";
        arrow.src = arrow.src.replace("right", "down");
        arrow.title = getLocalizedString('click_to_collapse');
    }
    else {
        // Currently expanded, so collapse it.
        block.style.display = "none";
        arrow.src = arrow.src.replace("down", "right");
        arrow.title = getLocalizedString('click_to_expand');
    }
    return false; // Make browser ignore href.
}

// ===================================================
// Create and uniquely name two levels of upward navigation buttons
// for Functions -- By Category pages

var top_button_count = 0;
var current_section_id = 0;

function addTopOfPageButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=\"#top_of_page\" onMouseOver=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_top_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_top_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_top_up.gif\"  alt=\"" +
        getLocalizedString('back_to_top_of_page') +
        "\" title=\"" +
        getLocalizedString('back_to_top_of_page') +
        "\" name=\"uparrow" +
        top_button_count +
        "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}


function updateSectionId(id) {
    current_section_id = id;
}


function addTopOfSectionButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=" +
        "\"#" + current_section_id + "\"" +
        " onMouseOver=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_section_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_section_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_section_up.gif\"  alt=\"" +
        getLocalizedString('back_to_top_of_section') +
        "\" title=\"" +
        getLocalizedString('back_to_top_of_section') +
        "\" name=\"uparrow" +
        top_button_count +
        "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}

// ===================================================
// Create and write to the document stream HTML for
// the link to the Doc Feedback Survey site.
//
// Doing this through a JavaScript function is necessary
// to work around the an issue with pages that are found
// through the search facility of the help browser--
//
// When found as the result of a search,
// the document that is displayed in the Help browser
// is actually a temporary document with a trivial URL
// such as "text://5", not an actual page location.
//
// But the Help browser inserts a <BASE> element at the beginning
// of each such temporary page, and the <BASE> element stores the
// actual location.
//
// So this function tests the URL of the document for the expression "text://"
// and if that expression is found, attempts to use the URL stored in
// the <BASE> element.

function writeDocFeedbackSurveyLink() {
    var queryexpression = document.location.href;

    if (queryexpression.search(/text:\/\//) != -1) {
        var baseelement = document.getElementsByTagName("BASE")[0];
        queryexpression = baseelement.href;
    }
    survey_url_yes = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-YES";
    survey_url_no = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-NO";

    code = '<div style="padding-right:10px" class="feedbackblock">' + '<strong>' + getLocalizedString('was_this_topic_helpful') + '</strong> <input type="button" value="' + getLocalizedString('yes') + '" onClick="openWindow(\'' + survey_url_yes + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '&nbsp;&nbsp;' + '<input type="button" value="' + getLocalizedString('no') + '" onClick="openWindow(\'' + survey_url_no + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '</div>';
    document.write(code);
}


// Utility function replacing openWindow function used by the web-site survey link code.
// In the help browser, the original code would create a blank window before loading the URL into the system browser.
function openWindow(url, width, height, options, name) {
    // ignore the arguments, except url
    document.location = url;
} // end function openWindow


// Utility function for linking to feedback survey, as of R2012b.
function openFeedbackWindow(url) {
    window.open(url, "_blank");
} // end function openFeedbackWindow


// ===================================================
// Workaround for G801125.
// This global object check tests for IE8 or lower.
if (document.all && !document.getElementsByClassName) {
    document.createElement("section");
}


// ===================================================
// Function reference pages

$(window).load(function () {
    // Perform breadcrumb check in window load, since all the images in the breadcrumb
    // need to be loaded for correct width calculations.
    //getBreadcrumb().mwBreadcrumb();
    //delay the expanding to enable the page to load completely.
    setTimeout(function () {
        expandCollapsedContent();
    }, 0);
    $('body').scrollspy({ target: '.offcanvas_nav', offset: getScrollTopAdjustment() + 1});
});

$(document).ready(function () {
    // this function is derived from an earlier version of jquery
    // and we are using it here for backwards compatability.
    $.browserInfo = function () {
        var ua = navigator.userAgent.toLowerCase();
        var match = /(chrome)[ \/]([\w.]+)/.exec(ua) ||
            /(webkit)[ \/]([\w.]+)/.exec(ua) ||
            /(opera)(?:.*version|)[ \/]([\w.]+)/.exec(ua) ||
            /(msie) ([\w.]+)/.exec(ua) ||
            ua.indexOf("compatible") < 0 && /(mozilla)(?:.*? rv:([\w.]+)|)/.exec(ua) ||
            [];

        var browserName = match ? match [1] : navigator.appName;
        var browserVersion = match ? match[2] : navigator.appVersion;
        return {
            name: browserName,
            version: browserVersion
        }
    };
	
	$('select').each(function () {
		var select = $(this);
		var selectedValue = select.find('option[selected]').val();

		if (selectedValue) {
		  select.val(selectedValue);
		} else {
		  select.prop('selectedIndex', 0);
		}
	});


    $(".nav_scrollspy").addClass("nav");
    $(".nav_scrollspy ul").addClass("nav");

    if ($.getParameterByName("browser") !== "F1help") {
        if ($.fn.setupToc) {
            $('div.toc_container_wrapper').setupToc();
        }
    }

    //Perform JS code which has any user visible impact first.
    //Check image sizes. Do not scale animated images or any images with hotspots.

    var searchCrumbContainer = $('#search_crumb_container');
    //old template
    if (searchCrumbContainer.length > 0) {
        $('#doc_center_content img:not(".animated-image, [usemap]"), #content_container2 img:not(".animated-image, [usemap]")').scaleImage();
    }

    $('#doc_center_content img.animated-image, #content_container2 img.animated-image').animateImage();


    $('.collapse').on('hidden.bs.collapse shown.bs.collapse', function () {
        var id = $(this).attr("id");
        var expandParent = $("[data-target='#" + id +"']");
        checkExpandAllLinkState(expandParent);
        triggerContentResize();
    });

    addSmoothScroll();

    $('#content_container .expandAllLink').click(function (e) {
        e.stopPropagation();
        var link = $(this);
        var allExpanders = link.closest('.expandableContent').find("[data-toggle='collapse']");
        if (link.data('allexpanded')) {
            doCollapse(allExpanders);
            setExpandAllLinkState(link, "collapsed");
        } else {
            doExpand(allExpanders);
            setExpandAllLinkState(link, "expanded");
        }
    });

    $('#expandAllPage').click(function () {
        var allExpanders = $("#content_container").find("[data-toggle='collapse']");
        if ($(this).data('allexpanded')) {
            doCollapseAllInPage(allExpanders);
            $(this).data('allexpanded', false);
            $(this).html(getLocalizedString('expand_all_in_page'));
        } else {
            doExpandAllInPage(allExpanders);
            $(this).data('allexpanded', true);
            $(this).html(getLocalizedString('collapse_all_in_page'));
        }
    });

    applySearchHighlight();

    $(window).bind('toc_resize', function () {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
    });

    $(window).bind('resize', function (e) {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
        $('.toc_container_wrapper').trigger('window_resize');
        refreshScrollSpy();
    });

    $(window).bind('content_resize', function (e) {
        $('.toc_container_wrapper').trigger('window_resize');
        refreshScrollSpy();
    });

    $(window).bind('intrnllnk_clicked', function (e) {
        expandCollapsedContent();
    });

    $('#search_crumb_container').width($('#content_container').width());
    getBreadcrumb().trigger('window_resize');

});

function setExpandAllLinkState(link, state) {
    if (state === 'expanded') {
        link.data('allexpanded', true);
        link.html(getLocalizedString('collapse_all'));
    } else if (state === 'collapsed') {
        link.data('allexpanded', false);
        link.html(getLocalizedString('expand_all'));
    }
}

function applySearchHighlight() {
    if ($.fn.highlight) {
        var highlighterCSSClass = ['highlight_01', 'highlight_02' , 'highlight_03', 'highlight_04', 'highlight_05'];
        var searchHighlightTerm = $.getParameterByName('searchHighlight');
        if (searchHighlightTerm.length > 0) {
            var searchHighlightArray = searchHighlightTerm.match(/"[^"]+"|\S+/g);
            $.each(searchHighlightArray, function (index, value) {
                var searchTerm = value.replace(/^"|"$/g, '');
                var cssClass = highlighterCSSClass[index % highlighterCSSClass.length];
                $("#doc_center_content").highlight(searchTerm,
                    {className: cssClass, wordsOnly: true});
                var elements = $("#doc_center_content").find("." + cssClass);

                setTimeout(function () {
                    expandHighlightedElements(elements);
                }, 50);
            });

            $(document).keyup(function (e) {
                if (e.which === 27) {
                    var classArray = $.map(highlighterCSSClass, function (value) {
                        return "." + value;
                    });
                    var highlightedEl = $(classArray.join(","));
                    $.each(highlightedEl, function () {
                        $(this).removeClass();
                    });
                }
            });
        }
    }
}

function expandHighlightedElements(elements) {
    $.each(elements, function () {
        if (!$(this).is(":visible")) {
            var collapsedParent = $(this).closest('.collapse');
            var id = collapsedParent.attr("id");
            var expandParent = $("[data-target='#" + id +"']");
            if (expandParent.length > 0) {
                prepareEltForExpansion(expandParent, true);
                doExpand(expandParent);
            }
        }
    });
}

//helper method to fetch the breadcrumb.
function getBreadcrumb() {
    var breadcrumb;
    if ($("#breadcrumbs").length != 0) {
        breadcrumb = $("#breadcrumbs");
    } else {
        breadcrumb = $(".breadcrumbs:first");
    }
    return breadcrumb;
}

function expandCollapsedContent() {
    if (location.hash.length > 0) {

        var target = getInternalLinkTarget(location.hash);
        if (target.length > 0) {
            var expandParent = getExpandParentForAnchorTarget(target);
            prepareEltForExpansion(expandParent);

            //scroll to the target first.
            var scrollParameter = getScrollParameter();
            var scrollTop = target.offset().top - getScrollTopAdjustment();
            $(scrollParameter).scrollTop(scrollTop);
            $('.anchor_hinting').removeClass('anchor_hinting');
            target.addClass('anchor_hinting');
            setTimeout(function () {
                target.removeClass('anchor_hinting');
            }, 5000);
            doExpand(expandParent);
        }
    }
}

function prepareEltForExpansion(elt, noAnimation) {
        doExpandNestedParent(elt, noAnimation);
    }

function getExpandParentForAnchorTarget(target) {
    var collapseToggle = target.parents("[data-toggle='collapse']");
    if (collapseToggle.length > 0) {
        return collapseToggle.first();
}
        return target;
    }

function addSmoothScroll() {
    $(".intrnllnk").each(function () {
        var hash = this.hash;
        var target = getInternalLinkTarget(hash);
        if (target.length > 0) {
            $(this).click(function (evt) {
                evt.preventDefault();
                var expandParent = getExpandParentForAnchorTarget(target);
                prepareEltForExpansion(expandParent);

                var scrollParameter = getScrollParameter();
                var scrollTop = target.offset().top - getScrollTopAdjustment();
                $('.anchor_hinting').removeClass('anchor_hinting');
                $(scrollParameter).animate({scrollTop: scrollTop}, 700, function () {
                    target.addClass('anchor_hinting');
                    setTimeout(function () {
                        target.removeClass('anchor_hinting');
                    }, 5000);
                    doExpand(expandParent);
                });
                if(history.pushState) {
                    history.pushState(null, null, hash);
                }
                else {
                    location.hash = hash;
                }
                refreshScrollSpy();
            })
        }
    });
}

function getInternalLinkTarget(hash) {
    //search for anchor with given hash as "name" atrribute value;
    var target = [];

    //Remove the first '#' character from the name attribute. Escape any special character from the name/id.
    var escapedHash = hash.substring(1).replace(/([;&,.+*~':"!^#$%@\[\]\(\)=>\|])/g, '\\$1');

    target = $("#" + escapedHash);
    return target;
}

// Assumes that the elt passed in has data-toggle=collapse, and
// data-target specifying the element which needs to be expanded/collapsed
function findExpandableContent(elt) {
    var target = elt.data('target');
    return $(target);
}

function triggerContentResize() {
    $(window).trigger('content_resize');
}

function doExpand(elt_array, noAnimation) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('show');
    });
}

function doCollapse(elt_array, noAnimation) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('hide');
    });
}

function checkExpandAllLinkState(elt) {
    //Check if the expandable elt is nested within another expandable elt.
    var expandableParent = elt.closest('.expandableContent');

    // If element is not nested, or there is not expand all link, return
    if (expandableParent.length === 0 || expandableParent.find('.expandAllLink').length === 0) {
        return;
    }
    var expandAllLink = expandableParent.find('.expandAllLink:first');
    var expandableChildren = expandableParent.find("[data-toggle='collapse']");

    if (elt.hasClass('collapsed')) {
        var allChildrenCollapsed = true;
        expandableChildren.each(function () {
            if (!$(this).hasClass("collapsed")) {
                allChildrenCollapsed = false;
            }
        });
        if (allChildrenCollapsed) {
            setExpandAllLinkState(expandAllLink, "collapsed");
        }
    } else {
        var allChildrenExpanded = true;
        expandableChildren.each(function () {
            if ($(this).hasClass("collapsed")) {
                allChildrenExpanded = false;
            }
        });

        if (allChildrenExpanded) {
            setExpandAllLinkState(expandAllLink, "expanded");
        }
    }
}

function doExpandNestedParent(elt, noAnimation) {
    var expandableParent = elt.parents("[data-toggle='collapse']").first();
    if (expandableParent.length > 0) {
        var expandable = findExpandableContent(expandableParent);
        expandable.collapse('show');
    }
}

function doExpandAllInPage(elt_array) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('show');
    });
    $.each($('#content_container .expandAllLink'), function (i, link) {
        setExpandAllLinkState($(link), "expanded");
    });
}

function doCollapseAllInPage(elt_array) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('hide');
    });
    $.each($('#content_container .expandAllLink'), function (i, link) {
        setExpandAllLinkState($(link), "collapsed");
    });
}

function getScrollParameter() {
    //On IE and FF, the slow scroll parameter is the HTML dom element. On webkit, it is the body.
    var browserInfo = $.browserInfo();
    var scrollParameter = (browserInfo.name === 'msie' || browserInfo.name === 'mozilla') ?
        'html' : 'body';
    return scrollParameter;
}

function getScrollTopAdjustment() {
    var scrollTop = 0;
    var searchCrumbContainer = $('#search_crumb_container');
    //old template
    if (searchCrumbContainer.length > 0) {
        if (searchCrumbContainer.css('position') === 'fixed') {
            scrollTop = 103;
        } else {
            scrollTop = 10;
        }
    } else {
        scrollTop = $('.sticky_header_container').height() + 10;
    }
    return scrollTop;
}

function refreshScrollSpy() {
    $('body#responsive_offcanvas').each(function () {
        $(this).data('bs.scrollspy').options.offset = getScrollTopAdjustment() + 1;
        var $spy = $(this).scrollspy('refresh');
    });
    $(window).scroll();
}

// Copyright 2002-2016 The MathWorks, Inc.
