# keel::bignum

Arbitrary-precision signed integer (`keel::BigInt`) and canonical exact rational (`keel::Rational`) arithmetic library built on [`keel::keel`](https://github.com/keelcpp/keel).

## Features

- **`keel::BigInt`** (`<keel/bigint.h>`):
  - Signed magnitude base-$2^{64}$ representation with 2-limb (128-bit) Small Buffer Optimization (SBO) inline on the stack.
  - Arithmetic (`+`, `-`, `*`, `/`, `%`, `<<`, `>>`), Knuth Algorithm D multi-limb division (`div_mod`), `gcd`, `lcm`, `pow_u64`, `sqrt_exact`.
  - Miller-Rabin primality testing (`is_prime`, deterministic below $2^{64}$) and Pollard-Brent rho integer factorization (`factorize`).
  - Correctly rounded IEEE-754 double conversion (`to_f64`), decimal parsing (`try_parse`), and decimal formatting (`to_string`).
- **`keel::Rational`** (`<keel/rational.h>`):
  - Always canonical (`den > 0`, `gcd(|num|, den) == 1`, zero is `0/1`).
  - Exact lossless conversion from finite IEEE-754 doubles (`from_f64_bits`) and scale-invariant correctly rounded conversion to double (`to_f64`).
  - Exact arithmetic, integer powers (`pow_u64`), parsing (`try_parse`), and formatting (`to_string`).
- Zero STL in public headers, compiles cleanly under `-fno-exceptions -fno-rtti`.

## Usage

```cpp
#include <keel/bignum.h>

keel::BigInt a = keel::BigInt::pow_u64(keel::BigInt(2), 128) - keel::BigInt(1);
auto factors = keel::BigInt::factorize(a);

keel::Rational r = keel::Rational(1, 3) + keel::Rational(1, 6); // 1/2
```
