// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { FetchHttpClient } from "./fetchHttpClient";
export class BrowserFetchHttpClient extends FetchHttpClient {
    prepareRequest(_httpRequest) {
        return Promise.resolve({});
    }
    processRequest(_operationResponse) {
        return Promise.resolve();
    }
    // eslint-disable-next-line @azure/azure-sdk/ts-apisurface-standardized-verbs
    fetch(input, init) {
        return fetch(input, init);
    }
}
//# sourceMappingURL=browserFetchHttpClient.js.map