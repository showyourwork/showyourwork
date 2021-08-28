// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { BaseRequestPolicy } from "./requestPolicy";
export const DefaultKeepAliveOptions = {
    enable: true
};
export function keepAlivePolicy(keepAliveOptions) {
    return {
        create: (nextPolicy, options) => {
            return new KeepAlivePolicy(nextPolicy, options, keepAliveOptions || DefaultKeepAliveOptions);
        }
    };
}
/**
 * KeepAlivePolicy is a policy used to control keep alive settings for every request.
 */
export class KeepAlivePolicy extends BaseRequestPolicy {
    /**
     * Creates an instance of KeepAlivePolicy.
     *
     * @param nextPolicy -
     * @param options -
     * @param keepAliveOptions -
     */
    constructor(nextPolicy, options, keepAliveOptions) {
        super(nextPolicy, options);
        this.keepAliveOptions = keepAliveOptions;
    }
    /**
     * Sends out request.
     *
     * @param request -
     * @returns
     */
    async sendRequest(request) {
        request.keepAlive = this.keepAliveOptions.enable;
        return this._nextPolicy.sendRequest(request);
    }
}
//# sourceMappingURL=keepAlivePolicy.js.map