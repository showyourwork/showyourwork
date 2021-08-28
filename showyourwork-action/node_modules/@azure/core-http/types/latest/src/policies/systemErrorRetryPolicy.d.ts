import { HttpOperationResponse } from "../httpOperationResponse";
import { WebResourceLike } from "../webResource";
import { BaseRequestPolicy, RequestPolicy, RequestPolicyFactory, RequestPolicyOptions } from "./requestPolicy";
export declare function systemErrorRetryPolicy(retryCount?: number, retryInterval?: number, minRetryInterval?: number, maxRetryInterval?: number): RequestPolicyFactory;
/**
 * @param retryCount - The client retry count.
 * @param retryInterval - The client retry interval, in milliseconds.
 * @param minRetryInterval - The minimum retry interval, in milliseconds.
 * @param maxRetryInterval - The maximum retry interval, in milliseconds.
 */
export declare class SystemErrorRetryPolicy extends BaseRequestPolicy {
    retryCount: number;
    retryInterval: number;
    minRetryInterval: number;
    maxRetryInterval: number;
    constructor(nextPolicy: RequestPolicy, options: RequestPolicyOptions, retryCount?: number, retryInterval?: number, minRetryInterval?: number, maxRetryInterval?: number);
    sendRequest(request: WebResourceLike): Promise<HttpOperationResponse>;
}
//# sourceMappingURL=systemErrorRetryPolicy.d.ts.map