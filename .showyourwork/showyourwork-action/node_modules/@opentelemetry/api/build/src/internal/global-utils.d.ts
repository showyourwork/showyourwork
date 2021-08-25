import type { DiagAPI } from '../api/diag';
import { ContextManager } from '../context/types';
import { DiagLogger } from '../diag';
import { TextMapPropagator } from '../propagation/TextMapPropagator';
import type { TracerProvider } from '../trace/tracer_provider';
export declare function registerGlobal<Type extends keyof OTelGlobalAPI>(type: Type, instance: OTelGlobalAPI[Type], diag: DiagAPI, allowOverride?: boolean): boolean;
export declare function getGlobal<Type extends keyof OTelGlobalAPI>(type: Type): OTelGlobalAPI[Type] | undefined;
export declare function unregisterGlobal(type: keyof OTelGlobalAPI, diag: DiagAPI): void;
declare type OTelGlobalAPI = {
    version: string;
    diag?: DiagLogger;
    trace?: TracerProvider;
    context?: ContextManager;
    propagation?: TextMapPropagator;
};
export {};
//# sourceMappingURL=global-utils.d.ts.map