# keel::bignum

Arbitrary-precision signed integer (`keel::bignum::BigInt`) and canonical exact rational (`keel::bignum::Rational`) arithmetic library built on [`keel::keel`](https://github.com/keelcpp/keel).

## Features

- **`keel::bignum::BigInt`** (`<keel/bigint.h>`):
  - Signed-magnitude base-$2^{64}$ representation with 2-limb (128-bit) Small Buffer Optimization (SBO) inline on the stack.
  - Arithmetic (`+`, `-`, `*` with Karatsuba above 2048 bits, `/`, `%`), alias-safe Knuth Algorithm D multi-limb division (`div_mod`), `gcd`, `lcm`, `pow_u64`.
  - Modular arithmetic & number theory: `mod_euclid`, `gcd_ext`, `inv_mod`, 4-bit windowed `pow_mod` (including negative exponents), `jacobi`, `isqrt`, `sqrt_exact`.
  - Infinite two's-complement bitwise operations (`~`, `&`, `|`, `^`, `<<`, `>>`, `test_bit`, `set_bit`, `clear_bit`, `flip_bit`, `trailing_zeros`, `popcount`) matching Python `int` semantics.
  - Base-2..36 & hexadecimal parsing/formatting (`try_parse`, `try_parse_radix`, `try_parse_hex`, `to_string`, `to_string_radix`, `to_hex`) and unsigned big-endian/little-endian binary serialization (`from_bytes_be`, `from_bytes_le`, `to_bytes_be`, `to_bytes_le`, `byte_length`).
  - Primality, factorization & CSPRNG generation: Miller-Rabin `is_prime` (deterministic below $2^{64}$), Pollard-Brent rho `factorize`, and OS-CSPRNG `random_bits`, `random_below`, `random_prime`.
  - Correctly rounded IEEE-754 double conversion (`to_f64`).
- **`keel::bignum::Rational`** (`<keel/rational.h>`):
  - Always canonical (`den > 0`, `gcd(|num|, den) == 1`, zero is `0/1`).
  - Exact lossless conversion from finite IEEE-754 doubles (`from_f64_bits`) and single-rounding scale-invariant conversion to double (`to_f64`, exact for normals and subnormals).
  - Exact arithmetic, `floor`, `ceil`, `trunc`, `inv`, `is_integer`, `pow_u64`, `pow_i64` (negative exponents), fraction parsing (`try_parse`), and decimal/scientific parsing (`try_parse_decimal`).
- Zero STL in public headers, compiles cleanly under `-fno-exceptions -fno-rtti`.

## Usage

```cpp
#include <keel/bignum.h>

using keel::bignum::BigInt;
using keel::bignum::Rational;

BigInt a = BigInt::pow_u64(BigInt(2), 128) - BigInt(1);
auto factors = BigInt::factorize(a);

auto p = BigInt::random_prime(256).value();
auto inv = BigInt::inv_mod(BigInt(65537), p - BigInt(1));
auto m = BigInt::pow_mod(BigInt(42), BigInt(65537), p).value();

Rational r = Rational(1, 3) + Rational(1, 6); // 1/2
auto r_neg = Rational::pow_i64(r, -3).value(); // 8
```

## Verification

- Unit test suite (`test_bignum`) + `keel` header purity/guidelines check (`test_guidelines_headers`).
- Differential verification against Python `int`, `math`, and `fractions.Fraction`:
  ```bash
  python tests/verify_reference.py
  ```
