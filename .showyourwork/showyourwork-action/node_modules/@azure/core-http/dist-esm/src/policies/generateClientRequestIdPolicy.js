// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { BaseRequestPolicy } from "./requestPolicy";
export function generateClientRequestIdPolicy(requestIdHeaderName = "x-ms-client-request-id") {
    return {
        create: (nextPolicy, options) => {
            return new GenerateClientRequestIdPolicy(nextPolicy, options, requestIdHeaderName);
        }
    };
}
export class GenerateClientRequestIdPolicy extends BaseRequestPolicy {
    constructor(nextPolicy, options, _requestIdHeaderName) {
        super(nextPolicy, options);
        this._requestIdHeaderName = _requestIdHeaderName;
    }
    sendRequest(request) {
        if (!request.headers.contains(this._requestIdHeaderName)) {
            request.headers.set(this._requestIdHeaderName, request.requestId);
        }
        return this._nextPolicy.sendRequest(request);
    }
}
//# sourceMappingURL=generateClientRequestIdPolicy.js.map