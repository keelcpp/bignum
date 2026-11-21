#include <keel/core.h>
#include <keel/math.h>
#include <keel/rational.h>

namespace keel {

Rational Rational::from_f64_bits(f64 v) {
    KEEL_ASSERT(keel::is_finite(v));
    i32 e = 0;
    const f64 f = keel::frexp(v, &e); // v == f * 2^e, f in [0.5, 1)
    // f has at most 53 significand bits (incl. hidden bit), so f * 2^53 is an exact integer.
    const i64 m = static_cast<i64>(f * 9007199254740992.0);
    if (m == 0) return Rational(); // +-0.0
    const i32 shift = e - 53;
    if (shift >= 0) {
        return Rational(BigInt(m) * BigInt::pow_u64(BigInt(i64{2}), static_cast<u64>(shift)), BigInt(i64{1}));
    }
    return Rational(BigInt(m), BigInt::pow_u64(BigInt(i64{2}), static_cast<u64>(-shift)));
}

f64 Rational::to_f64() const {
    if (m_num.is_zero()) return 0.0;
    if (m_den == BigInt(i64{1})) return m_num.to_f64();
    const usize bn = m_num.bit_length();
    const usize bd = m_den.bit_length();
    if (bn <= 53 && bd <= 53) {
        const f64 v = static_cast<f64>(m_num.low_limb()) / static_cast<f64>(m_den.low_limb());
        return m_num.is_negative() ? -v : v;
    }
    const BigInt n_abs = m_num.abs();
    const BigInt& d_val = m_den;
    i32 shift = 54 - (static_cast<i32>(bn) - static_cast<i32>(bd));
    BigInt ns = shift >= 0 ? (n_abs << static_cast<usize>(shift)) : n_abs;
    BigInt ds = shift < 0 ? (d_val << static_cast<usize>(-shift)) : d_val;
    if (ns < (ds << 54)) {
        ns = ns << 1;
        ++shift;
    } else if (ns >= (ds << 55)) {
        ds = ds << 1;
        --shift;
    }
    BigInt q;
    BigInt r;
    BigInt::div_mod(ns, ds, q, r);
    const u64 q64 = q.low_limb();
    u64 mant53 = q64 >> 2;
    const bool guard = ((q64 >> 1) & 1ULL) != 0;
    const bool sticky = ((q64 & 1ULL) != 0) || !r.is_zero();
    if (guard && (sticky || (mant53 & 1ULL) != 0)) ++mant53;
    const f64 scaled = keel::ldexp(static_cast<f64>(mant53), 2 - shift);
    return m_num.is_negative() ? -scaled : scaled;
}

Rational::Rational(i64 v) noexcept : m_num(v), m_den(i64{1}) {}

Rational::Rational(u64 v) noexcept : m_num(v), m_den(i64{1}) {}

Rational::Rational(BigInt p, BigInt q) : m_num(keel::move(p)), m_den(keel::move(q)) { canonicalize(); }

void Rational::canonicalize() {
    KEEL_ASSERT(!m_den.is_zero());
    if (m_num.is_zero()) {
        m_den = BigInt(i64{1});
        return;
    }
    if (m_den.is_negative()) {
        m_num = -m_num;
        m_den = -m_den;
    }
    const BigInt g = BigInt::gcd(m_num, m_den);
    if (g > BigInt(i64{1})) {
        m_num = m_num / g;
        m_den = m_den / g;
    }
}

Rational Rational::operator-() const {
    Rational out;
    out.m_num = -m_num;
    out.m_den = m_den;
    return out;
}

Rational Rational::abs() const {
    Rational out = *this;
    out.m_num = m_num.abs();
    return out;
}

Rational operator+(const Rational& a, const Rational& b) {
    return Rational(a.num() * b.den() + b.num() * a.den(), a.den() * b.den());
}

Rational operator-(const Rational& a, const Rational& b) {
    return Rational(a.num() * b.den() - b.num() * a.den(), a.den() * b.den());
}

Rational operator*(const Rational& a, const Rational& b) {
    return Rational(a.num() * b.num(), a.den() * b.den());
}

Rational operator/(const Rational& a, const Rational& b) {
    KEEL_ASSERT(!b.is_zero());
    return Rational(a.num() * b.den(), a.den() * b.num());
}

bool operator==(const Rational& a, const Rational& b) noexcept {
    return a.num() == b.num() && a.den() == b.den();
}

bool operator!=(const Rational& a, const Rational& b) noexcept { return !(a == b); }

bool operator<(const Rational& a, const Rational& b) { return a.num() * b.den() < b.num() * a.den(); }

bool operator<=(const Rational& a, const Rational& b) { return !(b < a); }

bool operator>(const Rational& a, const Rational& b) { return b < a; }

bool operator>=(const Rational& a, const Rational& b) { return !(a < b); }

Rational Rational::pow_u64(const Rational& base, u64 exp) {
    Rational result(i64{1});
    Rational b = base;
    while (exp != 0) {
        if (exp & 1ULL) result = result * b;
        exp >>= 1;
        if (exp != 0) b = b * b;
    }
    return result;
}

u64 Rational::hash_bits() const noexcept {
    return (num().hash_bits() * 0x9E3779B97F4A7C15ULL) ^ (den().hash_bits() + 0x165667B19E3779F9ULL);
}

bool Rational::fits_i64(i64& out) const noexcept {
    if (den() != BigInt(i64{1})) return false;
    return num().fits_i64(out);
}

String Rational::to_string(Allocator* alloc) const {
    if (den() == BigInt(i64{1})) return num().to_string(alloc);
    String out = num().to_string(alloc);
    out.push_back('/');
    out.append(den().to_string(alloc));
    return out;
}

bool Rational::try_parse(Rational& out, StringRef s) {
    out = Rational();
    const usize slash = s.find('/');
    if (slash == kInvalidIndex) {
        BigInt v;
        if (!BigInt::try_parse(v, s)) return false;
        out = Rational(keel::move(v), BigInt(i64{1}));
        return true;
    }
    BigInt p;
    BigInt q;
    if (!BigInt::try_parse(p, s.substr(0, slash))) return false;
    if (!BigInt::try_parse(q, s.substr(slash + 1))) return false;
    if (q.is_zero() || q.is_negative()) return false;
    out = Rational(keel::move(p), keel::move(q));
    return true;
}

} // namespace keel
