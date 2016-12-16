// Copyright 2014-2015 The MathWorks, Inc.
/**
    Wrapper for StandaloneEqnRenderer which takes care of loading dojo.

    Usage: just include this script from HTML. By default it selects all <math> elements
           on a page and replaces them by the rendered output.

           Optionally other configuration settings for StandaloneEqnRenderer can be passed
           through by setting a configuration object in the variable "equationrendererConfig"
           (like dojoConfig for dojo):

        <script>
            var equationrendererConfig = {
                equationEncoding: "attribute"
            }
        </script>
        <script type="text/javascript"
                src="<matlab root path>/ui/equationrenderer/release/MathRenderer.js"></script>

 **/

(function () {
    var isSupportedBrowser = function () {
        // IE 8 and prior doesn't have this function.
        return (document.addEventListener);
    };

    var replaceMathTagsByAltimgs = function () {
        // Use altimg attributes of math tags as fallback.
        var mathTags = document.getElementsByTagName("math");
        var idx, math, altimg, domImg, altimgWidth, altimgHeight, alttext, altimgValign, span;
        for (idx = mathTags.length - 1; idx >= 0; idx = idx - 1) {
            math = mathTags[idx];
            altimg = math.altimg || "";
            if (altimg) {
                domImg = document.createElement("img");
                domImg.src = altimg;
                altimgWidth = math["altimg-width"] || 0;
                if (altimgWidth) {
                    domImg.width = altimgWidth;
                }
                altimgHeight = math["altimg-height"] || 0;
                if (altimgHeight) {
                    domImg.height = altimgHeight;
                }
                alttext = math.alttext || "";
                if (alttext) {
                    domImg.alt = alttext;
                }

                altimgValign = math["altimg-valign"] || "";
                span = document.createElement("span");
                span.appendChild(domImg);
                if (altimgValign) {
                    span.style.cssText = "vertical-align: " + altimgValign + "; ";
                }

                // Replace math in DOM by the image tag.
                math.parentNode.replaceChild(span, math);
            } else {
                // Just remove the math tag.
                math.parentNode.removeChild(math);
            }
        }
    };

    if (!isSupportedBrowser()) {
        window.onload = replaceMathTagsByAltimgs;
        return;
    }

    // These are the configuration parameters which are passed through to the
    // StandaloneEqnRenderer by default.
    var equationrendererDefaultConfig = {
            flavor: "MathType",
            equationFormat: "mathml",
            equationEncoding: "element",
            equationRootElement: "math",
            cacheFontMetrics: false
        };

    // global value, otherwise dojo doesn't see it
    this.dojoConfig = {
        isDebug: true,
        async: true,
        cacheBust: false,
        packages: [
            {name: "dojo", location: "./"},
            {name: "MW", location: "./MW"}
        ]
    };

    (function () {
        var userConfig = this.equationrendererConfig || {};
        var key;
        for (key in userConfig) {
            if (userConfig.hasOwnProperty(key)) {
                equationrendererDefaultConfig[key] = userConfig[key];
            }
        }

        var head = document.getElementsByTagName("head")[0];
        var style;

        // first hide existing math elements to avoid flickering equations
        if (equationrendererDefaultConfig.equationEncoding === "element") {
            style = document.createElement("style");
            style.type = "text/css";
            style.textContent = equationrendererDefaultConfig.equationRootElement +
                " { visibility: hidden; }";
            head.appendChild(style);
        }

        // load dojo first
        var script = document.createElement("script");
        var callback = function () {
            require([
                "MW/equations/renderer/StandaloneEqnRenderer", "dojo/domReady!"
            ], function (StandaloneEqnRenderer) {
                var ren = new StandaloneEqnRenderer(equationrendererDefaultConfig);
                setTimeout(function () {
                    ren.render();
                }, 0);
            });
        };
        script.type = "text/javascript";
        script.onreadystatechange = function () {
            if (this.readyState === "complete") {
                callback();
            }
        };
        // just in case, if script.onreadystatechange is not called
        script.onload = callback;

        // find the script tag which loads this file and extract the path
        var basePath = "";
        var scriptTags = (document.documentElement || document).getElementsByTagName("script");
        var thisFileName = new RegExp("(^.*)enderer.js$");
        var i, srcStr;
        for (i = scriptTags.length - 1; i >= 0; i = i - 1) {
            srcStr = scriptTags[i].getAttribute("src") || "";
            if (srcStr.match(thisFileName)) {
                // pick up the begin of the path and use it below to load dojo
                basePath = srcStr.slice(0, -15);
                break;
            }
        }

        // This code can be removed, once we use xstyle to load CSS.
        var link = document.createElement("link");
        link.type = "text/css";
        link.rel = "stylesheet";
        link.href = basePath + "css/main.css";
        head.appendChild(link);

        script.src = basePath + "dojo/dojo.js";
        head.appendChild(script);
    }());
}());
