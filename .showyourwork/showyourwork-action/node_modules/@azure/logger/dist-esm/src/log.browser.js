// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { __read, __spread } from "tslib";
export function log() {
    var args = [];
    for (var _i = 0; _i < arguments.length; _i++) {
        args[_i] = arguments[_i];
    }
    if (args.length > 0) {
        var firstArg = String(args[0]);
        if (firstArg.includes(":error")) {
            console.error.apply(console, __spread(args));
        }
        else if (firstArg.includes(":warning")) {
            console.warn.apply(console, __spread(args));
        }
        else if (firstArg.includes(":info")) {
            console.info.apply(console, __spread(args));
        }
        else if (firstArg.includes(":verbose")) {
            console.debug.apply(console, __spread(args));
        }
        else {
            console.debug.apply(console, __spread(args));
        }
    }
}
//# sourceMappingURL=log.browser.js.map