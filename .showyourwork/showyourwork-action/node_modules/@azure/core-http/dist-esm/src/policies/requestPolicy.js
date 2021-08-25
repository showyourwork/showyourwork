// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
import { HttpPipelineLogLevel } from "../httpPipelineLogLevel";
export class BaseRequestPolicy {
    constructor(_nextPolicy, _options) {
        this._nextPolicy = _nextPolicy;
        this._options = _options;
    }
    /**
     * Get whether or not a log with the provided log level should be logged.
     * @param logLevel - The log level of the log that will be logged.
     * @returns Whether or not a log with the provided log level should be logged.
     */
    shouldLog(logLevel) {
        return this._options.shouldLog(logLevel);
    }
    /**
     * Attempt to log the provided message to the provided logger. If no logger was provided or if
     * the log level does not meat the logger's threshold, then nothing will be logged.
     * @param logLevel - The log level of this log.
     * @param message - The message of this log.
     */
    log(logLevel, message) {
        this._options.log(logLevel, message);
    }
}
/**
 * Optional properties that can be used when creating a RequestPolicy.
 */
export class RequestPolicyOptions {
    constructor(_logger) {
        this._logger = _logger;
    }
    /**
     * Get whether or not a log with the provided log level should be logged.
     * @param logLevel - The log level of the log that will be logged.
     * @returns Whether or not a log with the provided log level should be logged.
     */
    shouldLog(logLevel) {
        return (!!this._logger &&
            logLevel !== HttpPipelineLogLevel.OFF &&
            logLevel <= this._logger.minimumLogLevel);
    }
    /**
     * Attempt to log the provided message to the provided logger. If no logger was provided or if
     * the log level does not meet the logger's threshold, then nothing will be logged.
     * @param logLevel - The log level of this log.
     * @param message - The message of this log.
     */
    log(logLevel, message) {
        if (this._logger && this.shouldLog(logLevel)) {
            this._logger.log(logLevel, message);
        }
    }
}
//# sourceMappingURL=requestPolicy.js.map