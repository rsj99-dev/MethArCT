# -*- coding: utf-8 -*-
"""A small, dependency-free JSON Schema validator.

Only the keyword subset used by the MethArCT AI schemas is implemented so that
report validation works without installing ``jsonschema``.  When the optional
``jsonschema`` package *is* installed the tests cross-check both
implementations.

Supported keywords: ``$ref`` (local ``#/definitions/...`` only), ``type``,
``required``, ``properties``, ``additionalProperties``, ``items``, ``enum``,
``const``, ``minimum``, ``maximum``, ``minItems``, ``maxItems``, ``pattern``,
``anyOf``, ``oneOf``.
"""

import re
from typing import Any, Dict, List, Optional

_TYPE_MAP = {
    "object": dict,
    "array": list,
    "string": str,
    "boolean": bool,
    "number": (int, float),
    "integer": int,
    "null": type(None),
}


def _type_matches(value: Any, expected: str) -> bool:
    expected_type = _TYPE_MAP.get(expected)
    if expected_type is None:
        return True
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected == "boolean":
        return isinstance(value, bool)
    return isinstance(value, expected_type)


def _resolve_ref(root: Dict[str, Any], ref: str) -> Dict[str, Any]:
    if not ref.startswith("#/"):
        raise ValueError("Only local JSON Pointer references are supported: {0}".format(ref))
    node = root
    for token in ref[2:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        if not isinstance(node, dict) or token not in node:
            raise ValueError("Unresolvable $ref: {0}".format(ref))
        node = node[token]
    return node


def validate(instance: Any, schema: Dict[str, Any], path: str = "$") -> List[str]:
    """Validate *instance* against *schema*; returns a list of error strings.

    The schema is treated as its own root document, so local ``$ref`` entries
    pointing into its ``definitions`` block resolve correctly.
    """
    if not isinstance(schema, dict):
        return []
    return validate_document(instance, schema)


def validate_document(instance: Any, document: Dict[str, Any]) -> List[str]:
    """Validate *instance* against a full schema *document* (with ``$ref`` support).

    This is the entry point used by the report / job / capabilities schemas,
    which place shared definitions under ``definitions``.
    """

    def _walk(value: Any, schema: Dict[str, Any], path: str) -> List[str]:
        if not isinstance(schema, dict):
            return []
        if "$ref" in schema:
            return _walk(value, _resolve_ref(document, schema["$ref"]), path)

        errors = []
        if "type" in schema:
            types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
            if not any(_type_matches(value, item) for item in types):
                errors.append(
                    "{0}: expected type {1}, got {2}".format(
                        path, "/".join(types), type(value).__name__
                    )
                )
                return errors

        if "const" in schema and value != schema["const"]:
            errors.append("{0}: expected constant {1!r}".format(path, schema["const"]))
        if "enum" in schema and value not in schema["enum"]:
            errors.append("{0}: value {1!r} not in {2!r}".format(path, value, schema["enum"]))

        if isinstance(value, str) and schema.get("pattern"):
            if re.search(schema["pattern"], value) is None:
                errors.append("{0}: string does not match pattern {1!r}".format(path, schema["pattern"]))

        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if schema.get("minimum") is not None and value < schema["minimum"]:
                errors.append("{0}: {1} < minimum {2}".format(path, value, schema["minimum"]))
            if schema.get("maximum") is not None and value > schema["maximum"]:
                errors.append("{0}: {1} > maximum {2}".format(path, value, schema["maximum"]))

        if isinstance(value, dict):
            for key in schema.get("required", []):
                if key not in value:
                    errors.append("{0}: missing required property {1!r}".format(path, key))
            properties = schema.get("properties", {})
            for key, item in value.items():
                if key in properties:
                    errors.extend(_walk(item, properties[key], "{0}.{1}".format(path, key)))
                else:
                    additional = schema.get("additionalProperties", True)
                    if additional is False:
                        errors.append("{0}: unexpected property {1!r}".format(path, key))
                    elif isinstance(additional, dict):
                        errors.extend(_walk(item, additional, "{0}.{1}".format(path, key)))

        if isinstance(value, list):
            if schema.get("minItems") is not None and len(value) < schema["minItems"]:
                errors.append("{0}: needs at least {1} item(s)".format(path, schema["minItems"]))
            if schema.get("maxItems") is not None and len(value) > schema["maxItems"]:
                errors.append("{0}: allows at most {1} item(s)".format(path, schema["maxItems"]))
            items = schema.get("items")
            if isinstance(items, dict):
                for index, item in enumerate(value):
                    errors.extend(_walk(item, items, "{0}[{1}]".format(path, index)))

        if "anyOf" in schema:
            if not any(not _walk(value, option, path) for option in schema["anyOf"]):
                errors.append("{0}: does not match any allowed schema".format(path))
        if "oneOf" in schema:
            matches = sum(1 for option in schema["oneOf"] if not _walk(value, option, path))
            if matches != 1:
                errors.append("{0}: matched {1} schemas, expected exactly 1".format(path, matches))

        return errors

    return _walk(instance, document, "$")


def first_error(instance: Any, document: Dict[str, Any]) -> Optional[str]:
    """Return the first validation error, or ``None`` when valid."""
    errors = validate_document(instance, document)
    return errors[0] if errors else None
