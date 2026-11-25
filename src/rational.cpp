#include <keel/core.h>
#include <keel/math.h>
#include <keel/rational.h>

namespace keel::bignum {

namespace {

constexpr i64 kMaxDecimalExponent = 100000;

} // namespace

Rational Rational::from_f64_bits(f64 v) {
    KEEL_ASSERT(keel::is_finite(v));
    if (!keel::is_finite(v)) return Rational();
    i32 e = 0;
    const f64 f = keel::frexp(v, &e); // v == f * 2^e, f in [0.5, 1)
    // f has at most 53 significand bits (incl. hidden bit), so f * 2^53 is an exact integer.
    const i64 m = as<i64>(f * 9007199254740992.0);
    if (m == 0) return Rational(); // +-0.0
    const i32 shift = e - 53;
    if (shift >= 0) {
        return Rational(BigInt(m) << as<usize>(shift), BigInt(i64{1}));
    }
    return Rational(BigInt(m), BigInt(i64{1}) << as<usize>(-shift));
}

f64 Rational::to_f64() const {
    if (m_num.is_zero()) return 0.0;
    if (m_den == BigInt(i64{1})) return m_num.to_f64();

    const usize bn = m_num.bit_length();
    const usize bd = m_den.bit_length();
    if (bn <= 53 && bd <= 53) {
        const f64 v = as<f64>(m_num.low_limb()) / as<f64>(m_den.low_limb());
        return m_num.is_negative() ? -v : v;
    }

    // Fast overflow / underflow bounds before any multi-limb shifting.
    if (bn > bd + 1025) {
        return m_num.is_negative() ? -Limits<f64>::infinity() : Limits<f64>::infinity();
    }
    if (bd > bn + 1076) {
        return m_num.is_negative() ? -0.0 : 0.0;
    }

    const BigInt n_abs = m_num.abs();
    const BigInt& d_val = m_den;

    // Find exact floor(log2(n_abs / d_val)), which is either (bn - bd) or (bn - bd - 1).
    i32 e = as<i32>(bn) - as<i32>(bd);
    if (e >= 0) {
        if (n_abs < (d_val << as<usize>(e))) --e;
    } else {
        if ((n_abs << as<usize>(-e)) < d_val) --e;
    }

    if (e >= 1024) {
        return m_num.is_negative() ? -Limits<f64>::infinity() : Limits<f64>::infinity();
    }
    if (e < -1075) {
        return m_num.is_negative() ? -0.0 : 0.0;
    }

    // Single rounding directly to target ULP (normal: e - 52; subnormal: -1074).
    const i32 lsb_exp = (e >= -1022) ? (e - 52) : -1074;
    const BigInt ns = (lsb_exp < 0) ? (n_abs << as<usize>(-lsb_exp)) : n_abs;
    const BigInt ds = (lsb_exp > 0) ? (d_val << as<usize>(lsb_exp)) : d_val;

    BigInt q;
    BigInt r;
    BigInt::div_mod(ns, ds, q, r);
    u64 mant = q.low_limb();

    // Exact round-to-nearest-even by comparing 2*r against ds.
    const BigInt two_r = r << 1;
    if (two_r > ds || (two_r == ds && (mant & 1ULL) != 0ULL)) {
        ++mant;
    }

    const f64 scaled = keel::ldexp(as<f64>(mant), lsb_exp);
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

BigInt Rational::trunc() const {
    if (is_integer()) return m_num;
    return m_num / m_den;
}

BigInt Rational::floor() const {
    if (is_integer()) return m_num;
    BigInt q;
    BigInt r;
    BigInt::div_mod(m_num, m_den, q, r);
    if (r.is_negative()) {
        q -= BigInt(i64{1});
    }
    return q;
}

BigInt Rational::ceil() const {
    if (is_integer()) return m_num;
    BigInt q;
    BigInt r;
    BigInt::div_mod(m_num, m_den, q, r);
    if (!r.is_zero() && !r.is_negative()) {
        q += BigInt(i64{1});
    }
    return q;
}

Optional<Rational> Rational::inv() const {
    if (m_num.is_zero()) return kNullOpt;
    Rational out;
    if (m_num.is_negative()) {
        out.m_num = -m_den;
        out.m_den = -m_num;
    } else {
        out.m_num = m_den;
        out.m_den = m_num;
    }
    return out;
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
    if (exp == 0ULL) return Rational(i64{1});
    if (base.is_zero()) return Rational();
    // Since gcd(|num|, den) == 1, gcd(|num|^exp, den^exp) == 1 (already canonical).
    Rational out;
    out.m_num = BigInt::pow_u64(base.m_num, exp);
    out.m_den = BigInt::pow_u64(base.m_den, exp);
    return out;
}

Optional<Rational> Rational::pow_i64(const Rational& base, i64 exp) {
    if (exp == 0) return Rational(i64{1});
    if (exp > 0) return pow_u64(base, as<u64>(exp));
    if (base.is_zero()) return kNullOpt;
    Optional<Rational> recip = base.inv();
    if (!recip) return kNullOpt;
    const u64 uexp = ~as<u64>(exp) + 1ULL;
    return pow_u64(*recip, uexp);
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

bool Rational::try_parse_decimal(Rational& out, StringRef s) {
    out = Rational();
    if (s.empty()) return false;

    usize i = 0;
    bool negative = false;
    if (s[i] == '+' || s[i] == '-') {
        negative = (s[i] == '-');
        ++i;
    }
    if (i >= s.size()) return false;

    String digits;
    digits.reserve(s.size());
    bool saw_dot = false;
    usize frac_digits = 0;

    while (i < s.size()) {
        const char ch = s[i];
        if (ch >= '0' && ch <= '9') {
            digits.push_back(ch);
            if (saw_dot) ++frac_digits;
            ++i;
        } else if (ch == '.' && !saw_dot) {
            saw_dot = true;
            ++i;
        } else {
            break;
        }
    }
    if (digits.empty()) return false;

    i64 exp10 = 0;
    if (i < s.size()) {
        if (s[i] != 'e' && s[i] != 'E') return false;
        ++i;
        bool exp_neg = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
            exp_neg = (s[i] == '-');
            ++i;
        }
        if (i >= s.size()) return false;
        while (i < s.size()) {
            const char ch = s[i];
            if (ch < '0' || ch > '9') return false;
            exp10 = exp10 * 10 + as<i64>(ch - '0');
            if (exp10 > kMaxDecimalExponent) return false;
            ++i;
        }
        if (exp_neg) exp10 = -exp10;
    }

    BigInt mant;
    if (!BigInt::try_parse(mant, digits)) return false;
    if (mant.is_zero()) {
        out = Rational();
        return true;
    }

    if (frac_digits > as<usize>(kMaxDecimalExponent)) return false;
    const i64 total_exp = exp10 - as<i64>(frac_digits);
    if (total_exp > kMaxDecimalExponent || total_exp < -kMaxDecimalExponent) return false;

    if (negative) mant = -mant;
    if (total_exp >= 0) {
        out = Rational(mant * BigInt::pow_u64(BigInt(i64{10}), as<u64>(total_exp)), BigInt(i64{1}));
    } else {
        out = Rational(keel::move(mant), BigInt::pow_u64(BigInt(i64{10}), as<u64>(-total_exp)));
    }
    return true;
}

} // namespace keel::bignum
