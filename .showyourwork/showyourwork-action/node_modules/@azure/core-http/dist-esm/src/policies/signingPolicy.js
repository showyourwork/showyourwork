// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { BaseRequestPolicy } from "./requestPolicy";
export function signingPolicy(authenticationProvider) {
    return {
        create: (nextPolicy, options) => {
            return new SigningPolicy(nextPolicy, options, authenticationProvider);
        }
    };
}
export class SigningPolicy extends BaseRequestPolicy {
    constructor(nextPolicy, options, authenticationProvider) {
        super(nextPolicy, options);
        this.authenticationProvider = authenticationProvider;
    }
    signRequest(request) {
        return this.authenticationProvider.signRequest(request);
    }
    sendRequest(request) {
        return this.signRequest(request).then((nextRequest) => this._nextPolicy.sendRequest(nextRequest));
    }
}
//# sourceMappingURL=signingPolicy.js.map