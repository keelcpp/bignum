#pragma once

#include "keel/bigint.h"
#include "keel/core.h"
#include "keel/string.h"

namespace keel::bignum {

// ============================================================================
// Rational — exact rational number in canonical form.
//
// Representation: m_num / m_den with BigInt parts.
// Invariants: m_den > 0, gcd(|m_num|, m_den) == 1, zero is 0/1.
// All arithmetic funnels through the (p, q) constructor, which enforces the
// invariants, so every Rational is canonical.
// ============================================================================
class Rational {
public:
    Rational() noexcept = default;
    Rational(i64 v) noexcept;
    Rational(u64 v) noexcept;

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsSigned<T> && !IsSame<RemoveCv<T>, i64>)
    Rational(T v) noexcept : Rational(as<i64>(v)) {}

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsUnsigned<T> && !IsSame<RemoveCv<T>, u64>)
    Rational(T v) noexcept : Rational(as<u64>(v)) {}

    // q must be non-zero (KEEL_ASSERT in debug). Moves sign to numerator and reduces by gcd.
    Rational(BigInt p, BigInt q);

    // Strict parse: "<bigint>" or "<signed bigint>/<unsigned bigint>".
    // Denominator zero or malformed input returns false.
    [[nodiscard]] static bool try_parse(Rational& out, StringRef s);

    // Exact lossless rational value of a finite IEEE-754 double (m * 2^e).
    [[nodiscard]] static Rational from_f64_bits(f64 v);

    // Correctly rounded IEEE-754 double (scale-invariant for arbitrary bit lengths).
    [[nodiscard]] f64 to_f64() const;

    [[nodiscard]] bool is_zero() const noexcept { return m_num.is_zero(); }
    [[nodiscard]] i32 signum() const noexcept { return m_num.signum(); }
    [[nodiscard]] const BigInt& num() const noexcept { return m_num; }
    [[nodiscard]] const BigInt& den() const noexcept { return m_den; }

    [[nodiscard]] bool fits_i64(i64& out) const noexcept; // denominator 1, numerator fits

    // Content hash mixing numerator and denominator.
    [[nodiscard]] u64 hash_bits() const noexcept;

    [[nodiscard]] static Rational pow_u64(const Rational& base, u64 exp);

    [[nodiscard]] String to_string(Allocator* alloc = nullptr) const; // "p" if den == 1, else "p/q"

    [[nodiscard]] Rational operator-() const;
    [[nodiscard]] Rational abs() const;

    friend Rational operator+(const Rational& a, const Rational& b);
    friend Rational operator-(const Rational& a, const Rational& b);
    friend Rational operator*(const Rational& a, const Rational& b);
    friend Rational operator/(const Rational& a, const Rational& b); // b != 0

    Rational& operator+=(const Rational& o) { return *this = *this + o; }
    Rational& operator-=(const Rational& o) { return *this = *this - o; }
    Rational& operator*=(const Rational& o) { return *this = *this * o; }
    Rational& operator/=(const Rational& o) { return *this = *this / o; }

    friend bool operator==(const Rational& a, const Rational& b) noexcept;
    friend bool operator!=(const Rational& a, const Rational& b) noexcept;
    friend bool operator<(const Rational& a, const Rational& b);
    friend bool operator<=(const Rational& a, const Rational& b);
    friend bool operator>(const Rational& a, const Rational& b);
    friend bool operator>=(const Rational& a, const Rational& b);

private:
    void canonicalize();

    BigInt m_num;
    BigInt m_den = BigInt(i64{1});
};

} // namespace keel::bignum
