// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import * as tough from "tough-cookie";
import * as http from "http";
import * as https from "https";
import node_fetch from "node-fetch";
import { FetchHttpClient } from "./fetchHttpClient";
import { createProxyAgent, isUrlHttps } from "./proxyAgent";
function getCachedAgent(isHttps, agentCache) {
    return isHttps ? agentCache.httpsAgent : agentCache.httpAgent;
}
export class NodeFetchHttpClient extends FetchHttpClient {
    constructor() {
        super(...arguments);
        this.proxyAgents = {};
        this.keepAliveAgents = {};
        this.cookieJar = new tough.CookieJar(undefined, { looseMode: true });
    }
    getOrCreateAgent(httpRequest) {
        const isHttps = isUrlHttps(httpRequest.url);
        // At the moment, proxy settings and keepAlive are mutually
        // exclusive because the 'tunnel' library currently lacks the
        // ability to create a proxy with keepAlive turned on.
        if (httpRequest.proxySettings) {
            let agent = getCachedAgent(isHttps, this.proxyAgents);
            if (agent) {
                return agent;
            }
            const tunnel = createProxyAgent(httpRequest.url, httpRequest.proxySettings, httpRequest.headers);
            agent = tunnel.agent;
            if (tunnel.isHttps) {
                this.proxyAgents.httpsAgent = tunnel.agent;
            }
            else {
                this.proxyAgents.httpAgent = tunnel.agent;
            }
            return agent;
        }
        else if (httpRequest.keepAlive) {
            let agent = getCachedAgent(isHttps, this.keepAliveAgents);
            if (agent) {
                return agent;
            }
            const agentOptions = {
                keepAlive: httpRequest.keepAlive
            };
            if (isHttps) {
                agent = this.keepAliveAgents.httpsAgent = new https.Agent(agentOptions);
            }
            else {
                agent = this.keepAliveAgents.httpAgent = new http.Agent(agentOptions);
            }
            return agent;
        }
        else {
            return isHttps ? https.globalAgent : http.globalAgent;
        }
    }
    // eslint-disable-next-line @azure/azure-sdk/ts-apisurface-standardized-verbs
    async fetch(input, init) {
        return node_fetch(input, init);
    }
    async prepareRequest(httpRequest) {
        const requestInit = {};
        if (this.cookieJar && !httpRequest.headers.get("Cookie")) {
            const cookieString = await new Promise((resolve, reject) => {
                this.cookieJar.getCookieString(httpRequest.url, (err, cookie) => {
                    if (err) {
                        reject(err);
                    }
                    else {
                        resolve(cookie);
                    }
                });
            });
            httpRequest.headers.set("Cookie", cookieString);
        }
        // Set the http(s) agent
        requestInit.agent = this.getOrCreateAgent(httpRequest);
        requestInit.compress = httpRequest.decompressResponse;
        return requestInit;
    }
    async processRequest(operationResponse) {
        if (this.cookieJar) {
            const setCookieHeader = operationResponse.headers.get("Set-Cookie");
            if (setCookieHeader !== undefined) {
                await new Promise((resolve, reject) => {
                    this.cookieJar.setCookie(setCookieHeader, operationResponse.request.url, { ignoreError: true }, (err) => {
                        if (err) {
                            reject(err);
                        }
                        else {
                            resolve();
                        }
                    });
                });
            }
        }
    }
}
//# sourceMappingURL=nodeFetchHttpClient.js.map