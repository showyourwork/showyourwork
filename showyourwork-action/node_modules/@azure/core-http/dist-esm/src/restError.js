// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { custom } from "./util/inspect";
import { Sanitizer } from "./util/sanitizer";
const errorSanitizer = new Sanitizer();
export class RestError extends Error {
    constructor(message, code, statusCode, request, response) {
        super(message);
        this.name = "RestError";
        this.code = code;
        this.statusCode = statusCode;
        this.request = request;
        this.response = response;
        Object.setPrototypeOf(this, RestError.prototype);
    }
    /**
     * Logging method for util.inspect in Node
     */
    [custom]() {
        return `RestError: ${this.message} \n ${errorSanitizer.sanitize(this)}`;
    }
}
RestError.REQUEST_SEND_ERROR = "REQUEST_SEND_ERROR";
RestError.PARSE_ERROR = "PARSE_ERROR";
//# sourceMappingURL=restError.js.map