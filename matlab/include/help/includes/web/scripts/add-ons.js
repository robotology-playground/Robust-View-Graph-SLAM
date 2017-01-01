(function ($) {
    $(document).ready(function () {
                      
        $.getParameterByName = function(name) {
            name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
            var regexS = "[\\?&]" + name + "=([^&#]*)";
            var regex = new RegExp(regexS, 'g');
            var results;
            var value = "";
            while (true) {
            results = regex.exec(window.location.href);
            if (results == null) {
                break;
            }
                value = decodeURIComponent(results[1].replace(/\+/g, " "));
            }
            return value;
        };

                      
        if ($.getParameterByName("client") === "addons") {
            var channel = new MW.postmessagechannel.PostMessageChannel(window, '__linkclk__');
            channel.connect(window.parent);

            function sendToParent(messageType, message, callback) {
                channel.send(messageType, message, callback);
            }


            setTimeout(function () {
                addEventHandlers();
            }, 0);

            function addEventHandlers() {
                $("#content_container a").on('click', function (event) {
                    event.preventDefault();
                    var href = $(this).get(0).href;
                    sendToParent("openExternalLink", { "url": href });

                });

                $('a[href^="matlab:"]').off('click');
                $('.intrnllnk').off('click');


                $('a[href^="matlab:"]').on('click', function (event) {
                    event.preventDefault();
                    var href = $(this).attr('href');
                    sendToParent("openMatlabLink", { "url": href });
                });

                $(".intrnllnk").each(function () {
                    var hash = this.hash;
                    var target = getInternalLinkTarget(hash);
                    if (target.length > 0) {
                        $(this).click(function (evt) {
                            evt.preventDefault();
                          
                            var expandParent = getExpandParentForAnchorTarget(target);
                            prepareEltForExpansion(expandParent);

                            var scrollTop = target.offset().top - 10;
                            $('.anchor_hinting').removeClass('anchor_hinting');
                            sendToParent("scrollTo",
                                { "scrollOffset": scrollTop },
                                { "onSuccess": function () {
                                    target.addClass('anchor_hinting');
                                    setTimeout(function () {
                                    target.removeClass('anchor_hinting');
                                }, 5000);
                                doExpand(expandParent);

                                } });
                        });
                    }
                });

            }

        }
    });
})(jQuery);

