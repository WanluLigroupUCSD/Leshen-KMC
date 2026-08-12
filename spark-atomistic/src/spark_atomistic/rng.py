# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 Philox state, derivation, counter, and midpoint52 mapping."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import struct
from typing import Any

from .errors import DomainFailure


MASK32 = 0xFFFFFFFF
MASK128 = (1 << 128) - 1
M0 = 0xD2511F53
M1 = 0xCD9E8D57
W0 = 0x9E3779B9
W1 = 0xBB67AE85
ALGORITHM = "Philox4x32-10:errata-1-midpoint52"


def _mul_hi_lo(a: int, b: int) -> tuple[int, int]:
    product = a * b
    return (product >> 32) & MASK32, product & MASK32


def philox4x32_10(counter: tuple[int, int, int, int],
                   key: tuple[int, int]) -> tuple[int, int, int, int]:
    c0, c1, c2, c3 = counter
    k0, k1 = key
    for _ in range(10):
        hi0, lo0 = _mul_hi_lo(M0, c0)
        hi1, lo1 = _mul_hi_lo(M1, c2)
        c0, c1, c2, c3 = ((hi1 ^ c1 ^ k0) & MASK32, lo1,
                          (hi0 ^ c3 ^ k1) & MASK32, lo0)
        k0 = (k0 + W0) & MASK32
        k1 = (k1 + W1) & MASK32
    return c0, c1, c2, c3


def words_to_uniform(a: int, b: int) -> float:
    if not (0 <= a <= MASK32 and 0 <= b <= MASK32):
        raise DomainFailure("INVALID_STATE", "state invalid", component="rng",
                            requirement="E1-DET-002-B")
    q = (a << 20) | (b >> 12)
    return (2 * q + 1) * (2.0 ** -53)


def _words_to_int(words: tuple[int, int, int, int] | list[int]) -> int:
    return sum(word << (32 * index) for index, word in enumerate(words))


def _int_to_words(value: int) -> tuple[int, int, int, int]:
    return tuple((value >> (32 * index)) & MASK32 for index in range(4))  # type: ignore[return-value]


@dataclass(slots=True)
class PhiloxStream:
    key: tuple[int, int]
    initial_counter: int = 0
    next_counter: int = 0
    consumed_blocks: int = 0
    consumed_uniforms: int = 0
    buffered_block: tuple[int, int, int, int] | None = None
    next_pair: int = 0

    def clone(self) -> "PhiloxStream":
        return PhiloxStream(self.key, self.initial_counter, self.next_counter,
                            self.consumed_blocks, self.consumed_uniforms,
                            self.buffered_block, self.next_pair)

    def commit_from(self, clone: "PhiloxStream") -> None:
        if clone.key != self.key or clone.initial_counter != self.initial_counter:
            raise DomainFailure("INVALID_STATE", "state invalid", component="rng",
                                requirement="E2-KMC-003")
        self.next_counter = clone.next_counter
        self.consumed_blocks = clone.consumed_blocks
        self.consumed_uniforms = clone.consumed_uniforms
        self.buffered_block = clone.buffered_block
        self.next_pair = clone.next_pair

    def uniform(self) -> float:
        if self.buffered_block is None:
            self.buffered_block = philox4x32_10(_int_to_words(self.next_counter), self.key)
            self.next_counter = (self.next_counter + 1) & MASK128
            self.consumed_blocks += 1
            self.next_pair = 0
        offset = 2 * self.next_pair
        result = words_to_uniform(self.buffered_block[offset], self.buffered_block[offset + 1])
        self.consumed_uniforms += 1
        if self.next_pair == 0:
            self.next_pair = 1
        else:
            self.buffered_block = None
            self.next_pair = 0
        return result

    def checkpoint(self) -> dict[str, Any]:
        return {"algorithm": ALGORITHM,
                "buffered_block": None if self.buffered_block is None else list(self.buffered_block),
                "consumed_blocks": self.consumed_blocks,
                "consumed_uniforms": self.consumed_uniforms,
                "initial_counter": list(_int_to_words(self.initial_counter)),
                "key": list(self.key),
                "next_counter": list(_int_to_words(self.next_counter)),
                "next_pair": self.next_pair}

    @classmethod
    def restore(cls, value: Any) -> "PhiloxStream":
        required = {"algorithm", "buffered_block", "consumed_blocks", "consumed_uniforms",
                    "initial_counter", "key", "next_counter", "next_pair"}
        if not isinstance(value, dict) or set(value) != required:
            raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                requirement="E2-RNG-002")
        key = value["key"]
        initial = value["initial_counter"]
        following = value["next_counter"]
        buffered = value["buffered_block"]
        arrays = [key, initial, following] + ([] if buffered is None else [buffered])
        if (value["algorithm"] != ALGORITHM or not isinstance(key, list) or len(key) != 2
                or any(not isinstance(item, list) or len(item) != 4 for item in arrays[1:])
                or any(type(word) is not int or not 0 <= word <= MASK32
                       for item in arrays for word in item)
                or type(value["consumed_blocks"]) is not int or value["consumed_blocks"] < 0
                or type(value["consumed_uniforms"]) is not int or value["consumed_uniforms"] < 0
                or value["next_pair"] not in {0, 1}):
            raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                requirement="E2-RNG-002")
        initial_integer = _words_to_int(initial)
        next_integer = _words_to_int(following)
        blocks = value["consumed_blocks"]
        uniforms = value["consumed_uniforms"]
        expected_buffered = uniforms % 2 == 1
        if (blocks != (uniforms + 1) // 2
                or next_integer != (initial_integer + blocks) & MASK128
                or (buffered is not None) != expected_buffered
                or value["next_pair"] != (1 if expected_buffered else 0)):
            raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                requirement="E2-RNG-002")
        if buffered is not None:
            previous = (next_integer - 1) & MASK128
            if tuple(buffered) != philox4x32_10(_int_to_words(previous), (key[0], key[1])):
                raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                    requirement="E2-RNG-002")
        return cls((key[0], key[1]), initial_integer, next_integer, blocks, uniforms,
                   None if buffered is None else tuple(buffered), value["next_pair"])


def _stream_from_digest(material: bytes) -> PhiloxStream:
    output = hashlib.sha256(material).digest()
    key = struct.unpack(">II", output[:8])
    counter_words = struct.unpack(">IIII", output[8:24])
    counter = _words_to_int(list(counter_words))
    return PhiloxStream(key, counter, counter)


def derive_trajectory_stream(seed: int) -> PhiloxStream:
    return _stream_from_digest(b"spark-trajectory-stream/2\0" + struct.pack(">Q", seed))


def derive_saddle_stream(seed: int, state_id: str, search_class: str,
                         search_index: int) -> PhiloxStream:
    state_bytes = state_id.encode("utf-8")
    class_bytes = search_class.encode("utf-8")
    material = (b"spark-saddle-substream/2\0" + struct.pack(">Q", seed)
                + struct.pack(">I", len(state_bytes)) + state_bytes
                + struct.pack(">I", len(class_bytes)) + class_bytes
                + struct.pack(">Q", search_index))
    return _stream_from_digest(material)


def derive_stream(seed: int, state_id: str, search_class: str,
                  search_index: int) -> PhiloxStream:
    if state_id == "trajectory" and search_class == "trajectory":
        return derive_trajectory_stream(seed)
    return derive_saddle_stream(seed, state_id, search_class, search_index)
