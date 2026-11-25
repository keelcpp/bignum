#!/usr/bin/env python3
"""
Differential verification suite for keel::bignum against Python 3's built-in
arbitrary-precision int, math, fractions.Fraction, decimal.Decimal, and IEEE-754 struct.

Usage:
    python tests/verify_reference.py [path/to/test_bignum.exe]
"""

import math
import os
import random
import struct
import subprocess
import sys
import tempfile
from decimal import Decimal
from fractions import Fraction
from pathlib import Path


def find_exe() -> Path:
    if len(sys.argv) > 1:
        p = Path(sys.argv[1])
        if p.exists():
            return p
    root = Path(__file__).resolve().parent.parent
    candidates = [
        root / "build" / "Release" / "test_bignum.exe",
        root / "build" / "Release" / "test_bignum",
        root / "build" / "Debug" / "test_bignum.exe",
        root / "build" / "test_bignum.exe",
        root / "build" / "test_bignum",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError("Could not locate test_bignum executable; pass path as argv[1].")


def c_divmod(a: int, b: int) -> tuple[int, int]:
    """C truncated division semantics: a == q*b + r, |r| < |b|, sign(r) == sign(a) or 0."""
    q = int(abs(a) // abs(b))
    if (a < 0) != (b < 0):
        q = -q
    r = a - q * b
    return q, r


def mod_euclid(a: int, m: int) -> int:
    am = abs(m)
    return a % am


def py_jacobi(a: int, n: int) -> int:
    assert n > 0 and n % 2 == 1
    if n == 1:
        return 1
    a %= n
    res = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def format_radix(n: int, base: int) -> str:
    digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    if n == 0:
        return "0"
    neg = n < 0
    m = -n if neg else n
    out = []
    while m > 0:
        m, r = divmod(m, base)
        out.append(digits[r])
    if neg:
        out.append("-")
    return "".join(reversed(out))


def f64_to_hex(x: float) -> str:
    return struct.pack(">d", x).hex()


def bigint_to_f64_hex(n: int) -> str:
    try:
        v = float(n)
    except OverflowError:
        v = float("-inf") if n < 0 else float("inf")
    return f64_to_hex(v)


def rat_to_f64_hex(r: Fraction) -> str:
    try:
        v = float(r)
    except OverflowError:
        v = float("-inf") if r < 0 else float("inf")
    # Preserve -0.0 sign on negative underflow
    if v == 0.0 and r < 0:
        v = -0.0
    return f64_to_hex(v)


def rat_str(r: Fraction) -> str:
    if r.denominator == 1:
        return str(r.numerator)
    return f"{r.numerator}/{r.denominator}"


def miller_rabin_py(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rand_signed(rng: random.Random, max_bits: int) -> int:
    bits = rng.randint(0, max_bits)
    v = rng.getrandbits(bits)
    return -v if rng.getrandbits(1) else v


def main() -> int:
    exe = find_exe()
    rng = random.Random(0xC0FFEE2026)

    commands: list[str] = []
    checks: list[tuple[str, object]] = []

    def add_exact(cmd: str, expected: str) -> None:
        commands.append(cmd)
        checks.append(("exact", expected))

    def add_custom(cmd: str, fn) -> None:
        commands.append(cmd)
        checks.append(("custom", fn))

    # 1. Edge-case + randomized arithmetic (+, -, *, div_mod, mod_euclid) across
    #    small, 64-bit, 128-bit SBO, 512-bit, and 2500-bit Karatsuba boundaries.
    edge_vals = [
        0, 1, -1, 2, -2, 3, -3,
        (1 << 63) - 1, -(1 << 63), (1 << 64) - 1, -(1 << 64) + 1, (1 << 64), -(1 << 64),
        (1 << 127) - 1, -(1 << 127), (1 << 128) - 1, (1 << 128), -(1 << 128),
        (1 << 192) - 1, (1 << 96) + 1,
    ]
    pairs = []
    for i in range(len(edge_vals)):
        pairs.append((edge_vals[i], edge_vals[(i * 7 + 3) % len(edge_vals)]))
    for bit_cap in (16, 64, 128, 256, 512, 1024, 2400, 3200):
        count = 20 if bit_cap <= 512 else 8
        for _ in range(count):
            pairs.append((rand_signed(rng, bit_cap), rand_signed(rng, bit_cap)))

    for a, b in pairs:
        add_exact(f"ADD {a} {b}", str(a + b))
        add_exact(f"SUB {a} {b}", str(a - b))
        add_exact(f"MUL {a} {b}", str(a * b))
        if b != 0:
            q, r = c_divmod(a, b)
            add_exact(f"DIVMOD {a} {b}", f"{q} {r}")
            add_exact(f"MODEUC {a} {b}", str(mod_euclid(a, b)))

    # 2. Shifts & two's-complement bitwise operations (&, |, ^, ~, bit info).
    for a, b in pairs[:80]:
        add_exact(f"AND {a} {b}", str(a & b))
        add_exact(f"OR {a} {b}", str(a | b))
        add_exact(f"XOR {a} {b}", str(a ^ b))
        add_exact(f"NOT {a}", str(~a))
        sh = rng.choice([0, 1, 7, 31, 32, 63, 64, 65, 127, 128, 200])
        add_exact(f"SHL {a} {sh}", str(a << sh))
        add_exact(f"SHR {a} {sh}", str(a >> sh))
        idx = rng.choice([0, 1, 2, 63, 64, 127, 128, 300])
        bl = abs(a).bit_length()
        byl = (bl + 7) // 8
        tz = (abs(a) & -abs(a)).bit_length() - 1 if a != 0 else 0
        pc = abs(a).bit_count()
        tb = (a >> idx) & 1
        add_exact(f"BITINFO {a} {idx}", f"{bl} {byl} {tz} {pc} {tb}")

    # 3. GCD, LCM, extended GCD, modular inverse, modular exponentiation, Jacobi, isqrt.
    for _ in range(60):
        bits = rng.choice([16, 64, 128, 256, 512])
        a = rand_signed(rng, bits)
        b = rand_signed(rng, bits)
        add_exact(f"GCD {a} {b}", str(math.gcd(a, b)))
        add_exact(f"LCM {a} {b}", str(math.lcm(a, b)))

        def check_gcdext(out: str, a_val=a, b_val=b) -> bool:
            parts = out.split()
            if len(parts) != 3:
                return False
            g, x, y = map(int, parts)
            return g == math.gcd(a_val, b_val) and a_val * x + b_val * y == g

        add_custom(f"GCDEXT {a} {b}", check_gcdext)

        m = rng.getrandbits(bits) + 1
        if math.gcd(a, m) == 1:
            add_exact(f"INVMOD {a} {m}", str(pow(a, -1, m)))
        else:
            add_exact(f"INVMOD {a} {m}", "NONE")

        # Positive exponent pow_mod (small, medium, and large RSA-like exponents).
        e = rng.getrandbits(rng.choice([0, 1, 8, 32, 64, 256]))
        add_exact(f"POWMOD {a} {e} {m}", str(pow(a, e, m)))

        # Negative exponent pow_mod.
        e_neg = -rng.randint(1, 50)
        if math.gcd(a, m) == 1:
            add_exact(f"POWMOD {a} {e_neg} {m}", str(pow(a, e_neg, m)))
        else:
            add_exact(f"POWMOD {a} {e_neg} {m}", "NONE")

        # Odd positive modulus for Jacobi symbol.
        n_odd = (rng.getrandbits(bits) | 1)
        add_exact(f"JACOBI {a} {n_odd}", str(py_jacobi(a, n_odd)))

        # isqrt (non-negative & negative).
        if a >= 0:
            add_exact(f"ISQRT {a}", str(math.isqrt(a)))
        else:
            add_exact(f"ISQRT {a}", "NONE")

    # 1024-bit RSA-scale pow_mod checks.
    for _ in range(5):
        m = rng.getrandbits(1024) | (1 << 1023) | 1
        a = rng.getrandbits(1024)
        e = rng.getrandbits(1024)
        add_exact(f"POWMOD {a} {e} {m}", str(pow(a, e, m)))

    # 4. Base 2..36 radix parsing & formatting.
    for base_in in (2, 3, 7, 8, 10, 16, 23, 32, 36):
        for base_out in (2, 5, 8, 10, 16, 32, 36):
            v = rand_signed(rng, rng.choice([0, 16, 64, 200, 500]))
            s_in = format_radix(v, base_in)
            add_exact(f"RADIX {s_in} {base_in} {base_out}", format_radix(v, base_out))

    # 5. Unsigned big-endian & little-endian binary serialization.
    for _ in range(50):
        nbytes = rng.choice([0, 1, 2, 7, 8, 9, 16, 32, 64, 128])
        raw = rng.randbytes(nbytes)
        min_len = rng.choice([0, 1, 8, 32, 64])
        hex_in = raw.hex() if raw else "-"

        v_be = int.from_bytes(raw, "big")
        out_len_be = max((v_be.bit_length() + 7) // 8, min_len)
        be_be = v_be.to_bytes(out_len_be, "big").hex() if out_len_be > 0 else "-"
        le_be = v_be.to_bytes(out_len_be, "little").hex() if out_len_be > 0 else "-"
        add_exact(f"BYTES_BE {hex_in} {min_len}", f"{v_be} {be_be} {le_be}")

        v_le = int.from_bytes(raw, "little")
        out_len_le = max((v_le.bit_length() + 7) // 8, min_len)
        be_le = v_le.to_bytes(out_len_le, "big").hex() if out_len_le > 0 else "-"
        le_le = v_le.to_bytes(out_len_le, "little").hex() if out_len_le > 0 else "-"
        add_exact(f"BYTES_LE {hex_in} {min_len}", f"{v_le} {be_le} {le_le}")

    # 6. Exact IEEE-754 double conversion for BigInt (including halfway ties & overflow).
    f64_bigints = [
        0, 1, -1, (1 << 53) - 1, (1 << 53), (1 << 53) + 1, (1 << 53) + 2, (1 << 53) + 3,
        (1 << 100) + (1 << 46), (1 << 100) + (1 << 46) + 1, (1 << 100) + (1 << 47) + (1 << 46),
        ((1 << 53) - 1) << 971,             # DBL_MAX
        (((1 << 53) - 1) << 971) + (1 << 970) - 1,  # just below halfway -> DBL_MAX
        (((1 << 53) - 1) << 971) + (1 << 970),      # exact halfway -> +inf
        1 << 1024, -(1 << 1024), 1 << 2000,
    ]
    for _ in range(40):
        f64_bigints.append(rand_signed(rng, rng.choice([54, 64, 128, 500, 1024])))
    for v in f64_bigints:
        add_exact(f"BIGF64 {v}", bigint_to_f64_hex(v))

    # 7. Rational arithmetic, pow_i64, floor/ceil/trunc, exact to_f64 (incl. subnormals!),
    #    from_f64_bits, and decimal parsing.
    for _ in range(60):
        na = rand_signed(rng, rng.choice([8, 64, 150]))
        da = rng.getrandbits(rng.choice([8, 64, 150])) + 1
        nb = rand_signed(rng, rng.choice([8, 64, 150]))
        db = rng.getrandbits(rng.choice([8, 64, 150])) + 1
        ra = Fraction(na, da)
        rb = Fraction(nb, db)
        sa, sb = rat_str(ra), rat_str(rb)
        add_exact(f"RAT_OP + {sa} {sb}", rat_str(ra + rb))
        add_exact(f"RAT_OP - {sa} {sb}", rat_str(ra - rb))
        add_exact(f"RAT_OP * {sa} {sb}", rat_str(ra * rb))
        if rb != 0:
            add_exact(f"RAT_OP / {sa} {sb}", rat_str(ra / rb))
        add_exact(f"RAT_OP < {sa} {sb}", "1" if ra < rb else "0")
        add_exact(f"RAT_ROUND {sa}", f"{math.floor(ra)} {math.ceil(ra)} {math.trunc(ra)}")
        add_exact(f"RAT_F64 {sa}", rat_to_f64_hex(ra))

        e = rng.randint(-6, 6)
        if ra == 0 and e < 0:
            add_exact(f"RAT_POW {sa} {e}", "NONE")
        else:
            add_exact(f"RAT_POW {sa} {e}", rat_str(ra ** e))

    # Subnormal & boundary Rational::to_f64 and Rational::from_f64_bits tests!
    special_rats = [
        Fraction(0, 1),
        Fraction(1, 1 << 1075),                    # exact half of DBL_TRUE_MIN -> 0.0
        Fraction((1 << 60) + 1, 1 << 1135),        # just above half -> DBL_TRUE_MIN
        Fraction(3, 1 << 1075),                    # 1.5 * DBL_TRUE_MIN -> ties to even 2 * DBL_TRUE_MIN
        Fraction(-3, 1 << 1075),                   # -1.5 * DBL_TRUE_MIN -> -2 * DBL_TRUE_MIN
        Fraction((1 << 52) - 1, 1 << 1074),        # max subnormal
        Fraction((1 << 53) - 1, 1 << 1075),        # halfway between max subnormal and DBL_MIN -> DBL_MIN
        Fraction(1, 1 << 1022),                    # DBL_MIN
        Fraction(((1 << 53) - 1) << 971, 1),       # DBL_MAX
        Fraction((1 << 1024), 1),                  # +inf
        Fraction(-(1 << 1024), 1),                 # -inf
        Fraction((1 << 2000) + 1, 3 * (1 << 2000) + 2),  # huge 2000-bit ratio ~ 1/3
        Fraction((1 << 100) + 1, (1 << 1150) + 3), # huge ratio landing deep in subnormal range!
    ]
    for r in special_rats:
        add_exact(f"RAT_F64 {rat_str(r)}", rat_to_f64_hex(r))

    special_f64_bits = [
        0x0000000000000000, 0x8000000000000000,
        0x0000000000000001, 0x8000000000000001,
        0x000FFFFFFFFFFFFF, 0x0010000000000000,
        0x3FF0000000000000, 0xBFF0000000000000,
        0x3FB999999999999A, 0x400921FB54442D18,
        0x7FEFFFFFFFFFFFFF, 0xFFEFFFFFFFFFFFFF,
    ]
    for b64 in special_f64_bits:
        val = struct.unpack(">d", struct.pack(">Q", b64))[0]
        add_exact(f"RAT_FROM_F64 {b64:016x}", rat_str(Fraction(val)))

    # Decimal parsing checks.
    dec_strings = [
        "0", "-0", "+0.0", "0.1", "-0.125", ".5", "5.", "123.456",
        "1e-10", "-2.5E+12", "9.99999999999999999999e-50", "+00042.00500e+3",
    ]
    for ds in dec_strings:
        add_exact(f"RAT_DEC {ds}", rat_str(Fraction(Decimal(ds))))

    # 8. Primality & CSPRNG random_prime verification (independent Python Miller-Rabin).
    for n_test in (0, 1, 2, 3, 4, 5, 17, 561, 1105, 1729, (1 << 61) - 1, (1 << 89) - 1, (1 << 127) - 1, (1 << 127) - 3):
        add_exact(f"PRIME {n_test}", "1" if miller_rabin_py(n_test) else "0")

    for bits in (2, 8, 16, 32, 64, 128, 256, 512):
        def check_prime(out: str, expected_bits=bits) -> bool:
            if out == "NONE":
                return False
            p = int(out)
            return p.bit_length() == expected_bits and miller_rabin_py(p)

        add_custom(f"GENPRIME {bits}", check_prime)

    # 9. Factorization differential checks (product == n, strictly ascending prime factors).
    fac_targets = [
        1, 2, 3, 4, 360, 561, 65537, (1 << 32) + 1, (1 << 64) - 1,
        1000000007 * 1000000009,
    ]
    for _ in range(15):
        fac_targets.append(rng.getrandbits(rng.choice([24, 48, 64, 80])) + 2)
    for n_val in fac_targets:
        def check_fac(out: str, target=n_val) -> bool:
            if target == 1:
                return out == "1"
            parts = out.split()
            prod = 1
            prev_p = 0
            for tok in parts:
                p_s, e_s = tok.split("^")
                p, e = int(p_s), int(e_s)
                if p <= prev_p or e <= 0 or not miller_rabin_py(p):
                    return False
                prev_p = p
                prod *= p ** e
            return prod == target

        add_custom(f"FACTORIZE {n_val}", check_fac)

    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".txt", delete=False) as tf:
        batch_path = Path(tf.name)
        tf.write("\n".join(commands) + "\n")

    try:
        proc = subprocess.run(
            [str(exe), "--batch", str(batch_path)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            check=False,
        )
    finally:
        batch_path.unlink(missing_ok=True)

    if proc.returncode != 0:
        print(f"ERROR: test_bignum --batch exited with {proc.returncode}")
        print(proc.stderr)
        return 1

    outputs = proc.stdout.splitlines()
    if len(outputs) != len(checks):
        print(f"ERROR: expected {len(checks)} output lines, got {len(outputs)}")
        return 1

    failures = 0
    for idx, (cmd, (kind, expected), actual) in enumerate(zip(commands, checks, outputs)):
        ok = (actual == expected) if kind == "exact" else bool(expected(actual))
        if not ok:
            failures += 1
            if failures <= 20:
                print(f"FAIL [{idx}] cmd={cmd}")
                print(f"  expected: {expected if kind == 'exact' else '<predicate>'}")
                print(f"  actual:   {actual}")

    if failures:
        print(f"verify_reference: {failures}/{len(checks)} FAILED")
        return 1

    print(f"verify_reference: {len(checks)}/{len(checks)} differential checks PASSED against Python {sys.version.split()[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
