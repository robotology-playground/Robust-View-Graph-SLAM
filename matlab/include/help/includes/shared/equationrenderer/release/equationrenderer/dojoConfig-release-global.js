var dojoConfig = {
    tlmSiblingOfDojo: false,
    isDebug: false,
    locale: "en-us",
    async: true,
    has: {
        "production": 1
    },
    setLocale: function (localesSupported, overrideLanguage) {
        var locale = (overrideLanguage || navigator.language || navigator.userLanguage).toLowerCase();
        var splitLocale = locale.split('-',2);
        var i = localesSupported.indexOf(locale);
        if (i >= 0) {
            dojoConfig.locale = localesSupported[i];
            return;
        }
        for(var splitIndex = 0;  splitIndex  < 2; splitIndex ++) {
            var testLocale = splitLocale[splitIndex] || splitLocale[0];
            if (testLocale === "zh" || testLocale === "cn") {
                break;
            }
            for(i = 0; i < localesSupported.length; i++) {
                if (localesSupported[i].split('-')[splitIndex] === testLocale) {
                    dojoConfig.locale = localesSupported[i];
                    return;
                }
            }
        }
        dojoConfig.locale = localesSupported[0];
    }
};

