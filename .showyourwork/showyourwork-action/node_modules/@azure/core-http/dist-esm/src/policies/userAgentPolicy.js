// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { HttpHeaders } from "../httpHeaders";
import { Constants } from "../util/constants";
import { getDefaultUserAgentKey, getPlatformSpecificData } from "./msRestUserAgentPolicy";
import { BaseRequestPolicy } from "./requestPolicy";
function getRuntimeInfo() {
    const msRestRuntime = {
        key: "core-http",
        value: Constants.coreHttpVersion
    };
    return [msRestRuntime];
}
function getUserAgentString(telemetryInfo, keySeparator = " ", valueSeparator = "/") {
    return telemetryInfo
        .map((info) => {
        const value = info.value ? `${valueSeparator}${info.value}` : "";
        return `${info.key}${value}`;
    })
        .join(keySeparator);
}
export const getDefaultUserAgentHeaderName = getDefaultUserAgentKey;
export function getDefaultUserAgentValue() {
    const runtimeInfo = getRuntimeInfo();
    const platformSpecificData = getPlatformSpecificData();
    const userAgent = getUserAgentString(runtimeInfo.concat(platformSpecificData));
    return userAgent;
}
export function userAgentPolicy(userAgentData) {
    const key = !userAgentData || userAgentData.key === undefined || userAgentData.key === null
        ? getDefaultUserAgentKey()
        : userAgentData.key;
    const value = !userAgentData || userAgentData.value === undefined || userAgentData.value === null
        ? getDefaultUserAgentValue()
        : userAgentData.value;
    return {
        create: (nextPolicy, options) => {
            return new UserAgentPolicy(nextPolicy, options, key, value);
        }
    };
}
export class UserAgentPolicy extends BaseRequestPolicy {
    constructor(_nextPolicy, _options, headerKey, headerValue) {
        super(_nextPolicy, _options);
        this._nextPolicy = _nextPolicy;
        this._options = _options;
        this.headerKey = headerKey;
        this.headerValue = headerValue;
    }
    sendRequest(request) {
        this.addUserAgentHeader(request);
        return this._nextPolicy.sendRequest(request);
    }
    addUserAgentHeader(request) {
        if (!request.headers) {
            request.headers = new HttpHeaders();
        }
        if (!request.headers.get(this.headerKey) && this.headerValue) {
            request.headers.set(this.headerKey, this.headerValue);
        }
    }
}
//# sourceMappingURL=userAgentPolicy.js.map