# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Canonical JSON, portable scalar domain, and immutable deep values."""

from __future__ import annotations

import hashlib
import json
import math
from types import MappingProxyType
from typing import Any, Iterable

from .errors import DomainFailure


JSON_INTEGER_LIMIT = (1 << 53) - 1


def _pairs_without_duplicates(pairs: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DomainFailure("INVALID_INPUT", "duplicate JSON object key",
                                component="json", requirement="E2-JSON-001")
        result[key] = value
    return result


def _reject_constant(token: str) -> None:
    raise DomainFailure("NONFINITE_RESULT", f"nonfinite JSON number: {token}",
                        component="json", requirement="E2-JSON-003")


def _valid_unicode_scalar_text(value: str) -> bool:
    return all(not (0xD800 <= ord(ch) <= 0xDFFF) for ch in value)


def validate_portable_value(value: Any, path: str = "$") -> None:
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, int):
        if not -JSON_INTEGER_LIMIT <= value <= JSON_INTEGER_LIMIT:
            raise DomainFailure("INVALID_INPUT", f"integer outside portable domain at {path}",
                                component="json", requirement="E2-JSON-002")
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise DomainFailure("NONFINITE_RESULT", f"nonfinite number at {path}",
                                component="json", requirement="E2-JSON-003")
        return
    if isinstance(value, str):
        if not _valid_unicode_scalar_text(value):
            raise DomainFailure("INVALID_INPUT", f"non-scalar Unicode at {path}",
                                component="json", requirement="E2-JSON-001")
        return
    if isinstance(value, list) or isinstance(value, tuple):
        for index, item in enumerate(value):
            validate_portable_value(item, f"{path}[{index}]")
        return
    if isinstance(value, dict) or isinstance(value, MappingProxyType):
        for key, item in value.items():
            if not isinstance(key, str) or not _valid_unicode_scalar_text(key):
                raise DomainFailure("INVALID_INPUT", f"invalid object key at {path}",
                                    component="json", requirement="E2-JSON-001")
            validate_portable_value(item, f"{path}.{key}")
        return
    raise DomainFailure("INVALID_INPUT", f"non-JSON value at {path}",
                        component="json", requirement="E2-JSON-001")


def parse_json(raw: str | bytes | bytearray) -> Any:
    try:
        if isinstance(raw, (bytes, bytearray)):
            encoded = bytes(raw)
            if encoded.startswith(b"\xef\xbb\xbf"):
                raise DomainFailure("INVALID_INPUT", "UTF-8 BOM is forbidden",
                                    component="json", requirement="E2-JSON-001")
            text = encoded.decode("utf-8", errors="strict")
        elif isinstance(raw, str):
            text = raw
            if text.startswith("\ufeff"):
                raise DomainFailure("INVALID_INPUT", "UTF-8 BOM is forbidden",
                                    component="json", requirement="E2-JSON-001")
        else:
            raise DomainFailure("INVALID_INPUT", "input must be UTF-8 text or bytes",
                                component="json", requirement="E2-JSON-001")
        value = json.loads(text, object_pairs_hook=_pairs_without_duplicates,
                           parse_constant=_reject_constant)
        validate_portable_value(value)
        return value
    except DomainFailure:
        raise
    except (UnicodeError, json.JSONDecodeError, RecursionError, MemoryError) as exc:
        raise DomainFailure("INVALID_INPUT", "malformed UTF-8 JSON",
                            component="json", requirement="E2-JSON-001") from exc


def _canonical_string(value: str) -> str:
    pieces = ['"']
    for character in value:
        codepoint = ord(character)
        if character == '"':
            pieces.append('\\"')
        elif character == "\\":
            pieces.append("\\\\")
        elif codepoint <= 0x1F:
            pieces.append(f"\\u{codepoint:04x}")
        else:
            pieces.append(character)
    pieces.append('"')
    return "".join(pieces)


def _canonical_float(value: float) -> str:
    if not math.isfinite(value):
        raise DomainFailure("NONFINITE_RESULT", "nonfinite value rejected",
                            component="json", requirement="E2-CAN-004")
    if value == 0.0:
        return "0"
    negative = value < 0.0
    source = repr(abs(value)).lower()
    if "e" not in source and source.endswith(".0"):
        source = source[:-2]
    if "e" in source:
        mantissa, exponent_text = source.split("e")
        exponent = int(exponent_text)
    else:
        mantissa, exponent = source, 0
    digits = mantissa.replace(".", "")
    decimal_position = (mantissa.find(".") if "." in mantissa else len(mantissa)) + exponent
    magnitude = abs(value)
    # E3-CAN-001: a real at or above 2^53 must remain syntactically distinct
    # from an integer token, otherwise the portable integer-domain parser
    # rejects bytes produced by this encoder itself.
    if 1e-6 <= magnitude < 2**53:
        if decimal_position <= 0:
            rendered = "0." + "0" * (-decimal_position) + digits
        elif decimal_position >= len(digits):
            rendered = digits + "0" * (decimal_position - len(digits))
        else:
            rendered = digits[:decimal_position] + "." + digits[decimal_position:]
    else:
        adjusted_exponent = decimal_position - 1
        rendered = digits[0]
        if len(digits) > 1:
            rendered += "." + digits[1:]
        rendered += "e" + ("+" if adjusted_exponent >= 0 else "") + str(adjusted_exponent)
    return ("-" if negative else "") + rendered


def _canonical_fragment(value: Any) -> str:
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if type(value) is int:
        return str(value)
    if isinstance(value, float):
        return _canonical_float(value)
    if isinstance(value, str):
        return _canonical_string(value)
    if isinstance(value, list) or isinstance(value, tuple):
        return "[" + ",".join(_canonical_fragment(item) for item in value) + "]"
    if isinstance(value, dict) or isinstance(value, MappingProxyType):
        return "{" + ",".join(
            _canonical_string(key) + ":" + _canonical_fragment(value[key])
            for key in sorted(value)
        ) + "}"
    raise DomainFailure("INVALID_INPUT", "value cannot be canonicalized",
                        component="json", requirement="E2-CAN-001")


def canonical_bytes(value: Any) -> bytes:
    normalized = deep_thaw(value)
    validate_portable_value(normalized)
    try:
        return _canonical_fragment(normalized).encode("utf-8")
    except (TypeError, ValueError, UnicodeError, RecursionError, MemoryError) as exc:
        raise DomainFailure("INVALID_INPUT", "value cannot be canonicalized",
                            component="json", requirement="E2-CAN-001") from exc


def canonical_text(value: Any) -> str:
    return canonical_bytes(value).decode("utf-8")


def digest(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def deep_freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: deep_freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(deep_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(deep_freeze(item) for item in value)
    return value


def deep_thaw(value: Any) -> Any:
    if isinstance(value, MappingProxyType) or isinstance(value, dict):
        return {key: deep_thaw(item) for key, item in value.items()}
    if isinstance(value, tuple) or isinstance(value, list):
        return [deep_thaw(item) for item in value]
    return value
