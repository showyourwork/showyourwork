// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
/**
 * Defines the default token refresh buffer duration.
 */
export const TokenRefreshBufferMs = 2 * 60 * 1000; // 2 Minutes
/**
 * Provides an {@link AccessTokenCache} implementation which clears
 * the cached {@link AccessToken}'s after the expiresOnTimestamp has
 * passed.
 *
 * @deprecated No longer used in the bearer authorization policy.
 */
export class ExpiringAccessTokenCache {
    /**
     * Constructs an instance of {@link ExpiringAccessTokenCache} with
     * an optional expiration buffer time.
     */
    constructor(tokenRefreshBufferMs = TokenRefreshBufferMs) {
        this.cachedToken = undefined;
        this.tokenRefreshBufferMs = tokenRefreshBufferMs;
    }
    setCachedToken(accessToken) {
        this.cachedToken = accessToken;
    }
    getCachedToken() {
        if (this.cachedToken &&
            Date.now() + this.tokenRefreshBufferMs >= this.cachedToken.expiresOnTimestamp) {
            this.cachedToken = undefined;
        }
        return this.cachedToken;
    }
}
//# sourceMappingURL=accessTokenCache.js.map