import { SerializerOptions } from "./util/serializer.common";
export declare class Serializer {
    readonly modelMappers: {
        [key: string]: any;
    };
    readonly isXML?: boolean | undefined;
    constructor(modelMappers?: {
        [key: string]: any;
    }, isXML?: boolean | undefined);
    validateConstraints(mapper: Mapper, value: unknown, objectName: string): void;
    /**
     * Serialize the given object based on its metadata defined in the mapper
     *
     * @param mapper - The mapper which defines the metadata of the serializable object
     * @param object - A valid Javascript object to be serialized
     * @param objectName - Name of the serialized object
     * @param options - additional options to deserialization
     * @returns A valid serialized Javascript object
     */
    serialize(mapper: Mapper, object: unknown, objectName?: string, options?: SerializerOptions): any;
    /**
     * Deserialize the given object based on its metadata defined in the mapper
     *
     * @param mapper - The mapper which defines the metadata of the serializable object
     * @param responseBody - A valid Javascript entity to be deserialized
     * @param objectName - Name of the deserialized object
     * @param options - Controls behavior of XML parser and builder.
     * @returns A valid deserialized Javascript object
     */
    deserialize(mapper: Mapper, responseBody: unknown, objectName: string, options?: SerializerOptions): any;
}
export interface MapperConstraints {
    InclusiveMaximum?: number;
    ExclusiveMaximum?: number;
    InclusiveMinimum?: number;
    ExclusiveMinimum?: number;
    MaxLength?: number;
    MinLength?: number;
    Pattern?: RegExp;
    MaxItems?: number;
    MinItems?: number;
    UniqueItems?: true;
    MultipleOf?: number;
}
export declare type MapperType = SimpleMapperType | CompositeMapperType | SequenceMapperType | DictionaryMapperType | EnumMapperType;
export interface SimpleMapperType {
    name: "Base64Url" | "Boolean" | "ByteArray" | "Date" | "DateTime" | "DateTimeRfc1123" | "Object" | "Stream" | "String" | "TimeSpan" | "UnixTime" | "Uuid" | "Number" | "any";
}
export interface CompositeMapperType {
    name: "Composite";
    className?: string;
    modelProperties?: {
        [propertyName: string]: Mapper;
    };
    additionalProperties?: Mapper;
    uberParent?: string;
    polymorphicDiscriminator?: PolymorphicDiscriminator;
}
export interface SequenceMapperType {
    name: "Sequence";
    element: Mapper;
}
export interface DictionaryMapperType {
    name: "Dictionary";
    value: Mapper;
}
export interface EnumMapperType {
    name: "Enum";
    allowedValues: any[];
}
export interface BaseMapper {
    /**
     * Name for the xml element
     */
    xmlName?: string;
    /**
     * Xml element namespace
     */
    xmlNamespace?: string;
    /**
     * Xml element namespace prefix
     */
    xmlNamespacePrefix?: string;
    /**
     * Determines if the current property should be serialized as an attribute of the parent xml element
     */
    xmlIsAttribute?: boolean;
    /**
     * Name for the xml elements when serializing an array
     */
    xmlElementName?: string;
    /**
     * Whether or not the current property should have a wrapping XML element
     */
    xmlIsWrapped?: boolean;
    /**
     * Whether or not the current property is readonly
     */
    readOnly?: boolean;
    /**
     * Whether or not the current property is a constant
     */
    isConstant?: boolean;
    /**
     * Whether or not the current property is required
     */
    required?: boolean;
    /**
     * Whether or not the current property allows mull as a value
     */
    nullable?: boolean;
    /**
     * The name to use when serializing
     */
    serializedName?: string;
    /**
     * Type of the mapper
     */
    type: MapperType;
    /**
     * Default value when one is not explicitly provided
     */
    defaultValue?: any;
    /**
     * Constraints to test the current value against
     */
    constraints?: MapperConstraints;
}
export declare type Mapper = BaseMapper | CompositeMapper | SequenceMapper | DictionaryMapper | EnumMapper;
export interface PolymorphicDiscriminator {
    serializedName: string;
    clientName: string;
    [key: string]: string;
}
export interface CompositeMapper extends BaseMapper {
    type: CompositeMapperType;
}
export interface SequenceMapper extends BaseMapper {
    type: SequenceMapperType;
}
export interface DictionaryMapper extends BaseMapper {
    type: DictionaryMapperType;
    headerCollectionPrefix?: string;
}
export interface EnumMapper extends BaseMapper {
    type: EnumMapperType;
}
export interface UrlParameterValue {
    value: string;
    skipUrlEncoding: boolean;
}
export declare function serializeObject(toSerialize: unknown): any;
export declare const MapperType: {
    Date: "Date";
    Base64Url: "Base64Url";
    Boolean: "Boolean";
    ByteArray: "ByteArray";
    DateTime: "DateTime";
    DateTimeRfc1123: "DateTimeRfc1123";
    Object: "Object";
    Stream: "Stream";
    String: "String";
    TimeSpan: "TimeSpan";
    UnixTime: "UnixTime";
    Number: "Number";
    Composite: "Composite";
    Sequence: "Sequence";
    Dictionary: "Dictionary";
    Enum: "Enum";
};
//# sourceMappingURL=serializer.d.ts.map