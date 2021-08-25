// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { __read, __spread } from "tslib";
import util from "util";
import { EOL } from "os";
export function log(message) {
    var args = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        args[_i - 1] = arguments[_i];
    }
    process.stderr.write("" + util.format.apply(util, __spread([message], args)) + EOL);
}
//# sourceMappingURL=log.js.map