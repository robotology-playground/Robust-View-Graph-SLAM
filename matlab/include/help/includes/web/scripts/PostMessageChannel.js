var MW;

if (MW === undefined) {
    MW = {};
}

if (MW.postmessagechannel === undefined) {
    MW.postmessagechannel = {};
}

MW.postmessagechannel.PostMessageChannel = (function () {
    "use strict";

    var Message, Identifier, PostMessageChannel;

    Identifier = {
        generate: function () {
            return new Date().getTime().toString();
        }
    };

    Message = function (type, body, identifier) {
        this.type = type;
        this.body = body;
        this.identifier = identifier || Identifier.generate();
    };

    Message.serialize = function (message, tag) {
        return tag + JSON.stringify([message.type, message.body, message.identifier]);
    };

    Message.deserialize = function (data, tag) {
        var typeBodyAndIdentifier;

        if (data.indexOf(tag) !== 0) {
            return new Message();
        }

        typeBodyAndIdentifier = JSON.parse(data.replace(tag, ""));

        return new Message(
            typeBodyAndIdentifier[0], typeBodyAndIdentifier[1], typeBodyAndIdentifier[2]
        );
    };

    Message.prototype.isValid = function () {
        return this.type && (this.body !== null);
    };

    Message.prototype.setEvent = function (event) {
        this.event = event;
    };

    /**
     * Lightweight cross-window communication over the postMessage API. Allows multiple window
     * objects to communicate (for example, a document containing an iframe can communicate
     * with the document in the iframe).
     */
    PostMessageChannel = function (localWindow, tag) {
        var self = this;

        // __mwpmc__ is hard-coded into R2015b MATLAB and should only be used by the add-ons desktop
        // integration layer.
        this.tag = tag || "__mwpmc__";

        this.localWindow = localWindow;

        this.listeners = {};

        this.callbacks = {};

        this.onMessage = function (event) {
            var message = Message.deserialize(event.data, self.tag);

            if (message.isValid() && self.hasListener(message.type)) {
                message.setEvent(event);
                self.listeners[message.type].call(self, message);
                if (message.type.match(/^__/) === null) {
                    self.send('__success', message.identifier);
                }
            }
        };

        this.on("__success", function (message) {
            var identifier = message.body;
            if(this.callbacks[identifier]) {
                if(this.callbacks[identifier].onSuccess) {
                    this.callbacks[identifier].onSuccess.call(this);
                    delete this.callbacks[identifier];
                }
            }
        });
        this.addPostMessageListener();
    };

    /**
     * Clears the channel's event listeners. Also removes channel-related listeners from
     * the local window object.
     */
    PostMessageChannel.prototype.disconnect = function () {
        this.targetWindow = null;
        this.listeners = {};
        this.localWindow.removeEventListener("message", this.onMessage);
    };

    /**
     * Waits for another channel to send a connect message. On receiving this message,
     * the message sender is set as the channel's target.
     */
    PostMessageChannel.prototype.listen = function () {
        this.on("__connect", function (message) {
            this.setTargetWindow(message.event.source);
        });
    };

    /**
     * Sends a connect message to another window, and sets that window as the channel's target.
     */
    PostMessageChannel.prototype.connect = function (targetWindow) {
        this.setTargetWindow(targetWindow);

        this.send("__connect", this.localWindow.location.origin);
        this.send("connect", this.localWindow.location.origin);
    };

    /**
     * Adds a callback for a specified message type. Within the callback :this: is bound
     * to the channel object.
     */
    PostMessageChannel.prototype.on = function (type, listener) {
        this.listeners[type] = listener;
    };

    PostMessageChannel.prototype.hasListener = function (type) {
        return this.listeners[type] !== undefined;
    };

    /**
     * Sends a message to the channel's target window.
     */
    PostMessageChannel.prototype.send = function (type, body, callbacks) {
        var message = new Message(type, body || "");
        if (callbacks !== undefined) {
            this.callbacks[message.identifier] = callbacks;
        }
        this.targetWindow.postMessage(Message.serialize(message, this.tag), "*");
    };

    PostMessageChannel.prototype.addPostMessageListener = function () {
        this.localWindow.addEventListener("message", this.onMessage, false);
    };

    PostMessageChannel.prototype.setTargetWindow = function (targetWindow) {
        this.targetWindow = targetWindow;
    };

    return PostMessageChannel;
}());
