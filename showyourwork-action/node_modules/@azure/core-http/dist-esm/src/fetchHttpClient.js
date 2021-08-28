// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { AbortController, AbortError } from "@azure/abort-controller";
import FormData from "form-data";
import { HttpHeaders } from "./httpHeaders";
import { RestError } from "./restError";
import { Transform } from "stream";
import { logger } from "./log";
export class ReportTransform extends Transform {
    constructor(progressCallback) {
        super();
        this.progressCallback = progressCallback;
        this.loadedBytes = 0;
    }
    _transform(chunk, _encoding, callback) {
        this.push(chunk);
        this.loadedBytes += chunk.length;
        this.progressCallback({ loadedBytes: this.loadedBytes });
        callback(undefined);
    }
}
export class FetchHttpClient {
    async sendRequest(httpRequest) {
        var _a;
        if (!httpRequest && typeof httpRequest !== "object") {
            throw new Error("'httpRequest' (WebResourceLike) cannot be null or undefined and must be of type object.");
        }
        const abortController = new AbortController();
        let abortListener;
        if (httpRequest.abortSignal) {
            if (httpRequest.abortSignal.aborted) {
                throw new AbortError("The operation was aborted.");
            }
            abortListener = (event) => {
                if (event.type === "abort") {
                    abortController.abort();
                }
            };
            httpRequest.abortSignal.addEventListener("abort", abortListener);
        }
        if (httpRequest.timeout) {
            setTimeout(() => {
                abortController.abort();
            }, httpRequest.timeout);
        }
        if (httpRequest.formData) {
            const formData = httpRequest.formData;
            const requestForm = new FormData();
            const appendFormValue = (key, value) => {
                // value function probably returns a stream so we can provide a fresh stream on each retry
                if (typeof value === "function") {
                    value = value();
                }
                if (value &&
                    Object.prototype.hasOwnProperty.call(value, "value") &&
                    Object.prototype.hasOwnProperty.call(value, "options")) {
                    requestForm.append(key, value.value, value.options);
                }
                else {
                    requestForm.append(key, value);
                }
            };
            for (const formKey of Object.keys(formData)) {
                const formValue = formData[formKey];
                if (Array.isArray(formValue)) {
                    for (let j = 0; j < formValue.length; j++) {
                        appendFormValue(formKey, formValue[j]);
                    }
                }
                else {
                    appendFormValue(formKey, formValue);
                }
            }
            httpRequest.body = requestForm;
            httpRequest.formData = undefined;
            const contentType = httpRequest.headers.get("Content-Type");
            if (contentType && contentType.indexOf("multipart/form-data") !== -1) {
                if (typeof requestForm.getBoundary === "function") {
                    httpRequest.headers.set("Content-Type", `multipart/form-data; boundary=${requestForm.getBoundary()}`);
                }
                else {
                    // browser will automatically apply a suitable content-type header
                    httpRequest.headers.remove("Content-Type");
                }
            }
        }
        let body = httpRequest.body
            ? typeof httpRequest.body === "function"
                ? httpRequest.body()
                : httpRequest.body
            : undefined;
        if (httpRequest.onUploadProgress && httpRequest.body) {
            const onUploadProgress = httpRequest.onUploadProgress;
            const uploadReportStream = new ReportTransform(onUploadProgress);
            if (isReadableStream(body)) {
                body.pipe(uploadReportStream);
            }
            else {
                uploadReportStream.end(body);
            }
            body = uploadReportStream;
        }
        const platformSpecificRequestInit = await this.prepareRequest(httpRequest);
        const requestInit = Object.assign({ body: body, headers: httpRequest.headers.rawHeaders(), method: httpRequest.method, signal: abortController.signal, redirect: "manual" }, platformSpecificRequestInit);
        let operationResponse;
        try {
            const response = await this.fetch(httpRequest.url, requestInit);
            const headers = parseHeaders(response.headers);
            const streaming = ((_a = httpRequest.streamResponseStatusCodes) === null || _a === void 0 ? void 0 : _a.has(response.status)) ||
                httpRequest.streamResponseBody;
            operationResponse = {
                headers: headers,
                request: httpRequest,
                status: response.status,
                readableStreamBody: streaming
                    ? response.body
                    : undefined,
                bodyAsText: !streaming ? await response.text() : undefined
            };
            const onDownloadProgress = httpRequest.onDownloadProgress;
            if (onDownloadProgress) {
                const responseBody = response.body || undefined;
                if (isReadableStream(responseBody)) {
                    const downloadReportStream = new ReportTransform(onDownloadProgress);
                    responseBody.pipe(downloadReportStream);
                    operationResponse.readableStreamBody = downloadReportStream;
                }
                else {
                    const length = parseInt(headers.get("Content-Length")) || undefined;
                    if (length) {
                        // Calling callback for non-stream response for consistency with browser
                        onDownloadProgress({ loadedBytes: length });
                    }
                }
            }
            await this.processRequest(operationResponse);
            return operationResponse;
        }
        catch (error) {
            const fetchError = error;
            if (fetchError.code === "ENOTFOUND") {
                throw new RestError(fetchError.message, RestError.REQUEST_SEND_ERROR, undefined, httpRequest);
            }
            else if (fetchError.type === "aborted") {
                throw new AbortError("The operation was aborted.");
            }
            throw fetchError;
        }
        finally {
            // clean up event listener
            if (httpRequest.abortSignal && abortListener) {
                let uploadStreamDone = Promise.resolve();
                if (isReadableStream(body)) {
                    uploadStreamDone = isStreamComplete(body);
                }
                let downloadStreamDone = Promise.resolve();
                if (isReadableStream(operationResponse === null || operationResponse === void 0 ? void 0 : operationResponse.readableStreamBody)) {
                    downloadStreamDone = isStreamComplete(operationResponse.readableStreamBody);
                }
                Promise.all([uploadStreamDone, downloadStreamDone])
                    .then(() => {
                    var _a;
                    (_a = httpRequest.abortSignal) === null || _a === void 0 ? void 0 : _a.removeEventListener("abort", abortListener);
                    return;
                })
                    .catch((e) => {
                    logger.warning("Error when cleaning up abortListener on httpRequest", e);
                });
            }
        }
    }
}
function isReadableStream(body) {
    return body && typeof body.pipe === "function";
}
function isStreamComplete(stream) {
    return new Promise((resolve) => {
        stream.on("close", resolve);
        stream.on("end", resolve);
        stream.on("error", resolve);
    });
}
export function parseHeaders(headers) {
    const httpHeaders = new HttpHeaders();
    headers.forEach((value, key) => {
        httpHeaders.set(key, value);
    });
    return httpHeaders;
}
//# sourceMappingURL=fetchHttpClient.js.map