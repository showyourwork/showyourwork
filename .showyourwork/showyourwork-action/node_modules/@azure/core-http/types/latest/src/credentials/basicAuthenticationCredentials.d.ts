import { WebResourceLike } from "../webResource";
import { ServiceClientCredentials } from "./serviceClientCredentials";
export declare class BasicAuthenticationCredentials implements ServiceClientCredentials {
    userName: string;
    password: string;
    authorizationScheme: string;
    /**
     * Creates a new BasicAuthenticationCredentials object.
     *
     * @param userName - User name.
     * @param password - Password.
     * @param authorizationScheme - The authorization scheme.
     */
    constructor(userName: string, password: string, authorizationScheme?: string);
    /**
     * Signs a request with the Authentication header.
     *
     * @param webResource - The WebResourceLike to be signed.
     * @returns The signed request object.
     */
    signRequest(webResource: WebResourceLike): Promise<WebResourceLike>;
}
//# sourceMappingURL=basicAuthenticationCredentials.d.ts.map