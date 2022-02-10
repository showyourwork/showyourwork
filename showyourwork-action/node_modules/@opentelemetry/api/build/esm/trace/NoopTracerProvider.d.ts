import { Tracer } from './tracer';
import { TracerProvider } from './tracer_provider';
/**
 * An implementation of the {@link TracerProvider} which returns an impotent
 * Tracer for all calls to `getTracer`.
 *
 * All operations are no-op.
 */
export declare class NoopTracerProvider implements TracerProvider {
    getTracer(_name?: string, _version?: string): Tracer;
}
//# sourceMappingURL=NoopTracerProvider.d.ts.map