#include <keel/algo.h>
#include <keel/bigint.h>
#include <keel/core.h>
#include <keel/math.h>
#include <keel/number.h>

namespace keel::bignum {

namespace {

constexpr u64 kPow10_19 = 10000000000000000000ULL; // 10^19 < 2^64

const u64 kPow10[20] = {
    1ULL,
    10ULL,
    100ULL,
    1000ULL,
    10000ULL,
    100000ULL,
    1000000ULL,
    10000000ULL,
    100000000ULL,
    1000000000ULL,
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL,
    100000000000000000ULL,
    1000000000000000000ULL,
    10000000000000000000ULL,
};

i32 nlz32(u32 x) noexcept {
    i32 n = 0;
    if (x == 0) return 32;
    if ((x >> 16) == 0) { n += 16; x <<= 16; }
    if ((x >> 24) == 0) { n += 8;  x <<= 8; }
    if ((x >> 28) == 0) { n += 4;  x <<= 4; }
    if ((x >> 30) == 0) { n += 2;  x <<= 2; }
    if ((x >> 31) == 0) { n += 1; }
    return n;
}

i32 nlz64(u64 x) noexcept {
    if (x == 0) return 64;
    const u32 hi = static_cast<u32>(x >> 32);
    return hi != 0 ? nlz32(hi) : 32 + nlz32(static_cast<u32>(x));
}

// Full 64x64 -> 128-bit product without intrinsics or __int128 (MSVC-portable).
void mulu64(u64 a, u64 b, u64* lo, u64* hi) noexcept {
    const u64 mask = 0xFFFFFFFFULL;
    const u64 a_lo = a & mask;
    const u64 a_hi = a >> 32;
    const u64 b_lo = b & mask;
    const u64 b_hi = b >> 32;
    const u64 ll = a_lo * b_lo;
    const u64 lh = a_lo * b_hi;
    const u64 hl = a_hi * b_lo;
    const u64 hh = a_hi * b_hi;
    const u64 mid = lh + hl;
    const u64 mid_carry = (mid < lh) ? 1ULL : 0ULL;
    const u64 lo_res = ll + (mid << 32);
    const u64 lo_carry = (lo_res < ll) ? 1ULL : 0ULL;
    *lo = lo_res;
    *hi = hh + (mid >> 32) + (mid_carry << 32) + lo_carry;
}

template <typename LimbContainer>
void to_digits32(const LimbContainer& limbs, Array<u32>& d) {
    d.clear();
    d.reserve(limbs.size() * 2);
    for (u64 l : limbs) {
        d.push_back(static_cast<u32>(l));
        d.push_back(static_cast<u32>(l >> 32));
    }
    while (!d.empty() && d.back() == 0) d.pop_back();
}

template <typename LimbContainer>
void from_digits32(const Array<u32>& d, LimbContainer& limbs) {
    limbs.assign((d.size() + 1) / 2, 0);
    for (usize i = 0; i < d.size(); ++i) {
        if (i & 1) {
            limbs[i >> 1] |= static_cast<u64>(d[i]) << 32;
        } else {
            limbs[i >> 1] |= d[i];
        }
    }
    while (!limbs.empty() && limbs.back() == 0) limbs.pop_back();
}

// Knuth algorithm D (Hacker's Delight divmnu) on base-2^32 digits.
// Digits are little-endian with no leading zero digit; u may be empty.
void div_mod_digits(const Array<u32>& u, const Array<u32>& v, Array<u32>& q, Array<u32>& r) {
    const usize m = u.size();
    const usize n = v.size();

    if (n > m) {
        q.clear();
        r = u.clone();
        return;
    }

    if (n == 1) {
        q.assign(m, u32{0});
        u64 k = 0;
        for (usize j = m; j-- > 0;) {
            const u64 t = (k << 32) | u[j];
            q[j] = static_cast<u32>(t / v[0]);
            k = t % v[0];
        }
        r.assign(1, static_cast<u32>(k));
        while (q.size() > 1 && q.back() == 0) q.pop_back();
        return;
    }

    const i32 s = nlz32(v[n - 1]);
    Array<u32> vn(n);
    Array<u32> un(m + 1, u32{0});

    for (usize i = n; i-- > 1;) {
        vn[i] = s ? (v[i] << s) | (v[i - 1] >> (32 - s)) : v[i];
    }
    vn[0] = v[0] << s;

    for (usize i = m; i-- > 1;) {
        un[i] = s ? (u[i] << s) | (u[i - 1] >> (32 - s)) : u[i];
    }
    un[0] = u[0] << s;
    un[m] = s ? u[m - 1] >> (32 - s) : 0;

    q.assign(m - n + 1, u32{0});

    for (usize j = m - n + 1; j-- > 0;) {
        const u64 num = (static_cast<u64>(un[j + n]) << 32) | un[j + n - 1];
        u64 qhat = num / vn[n - 1];
        u64 rhat = num % vn[n - 1];

        while (qhat >= (1ULL << 32) || qhat * vn[n - 2] > ((rhat << 32) | un[j + n - 2])) {
            --qhat;
            rhat += vn[n - 1];
            if (rhat >= (1ULL << 32)) break;
        }

        u64 k = 0;
        for (usize i = 0; i < n; ++i) {
            const u64 p = qhat * vn[i];
            const u64 sub = k + (p & 0xFFFFFFFFULL);
            const u64 cur = un[i + j];
            un[i + j] = static_cast<u32>(cur - sub);
            const u64 deficit = sub > cur ? sub - cur : 0;
            k = (p >> 32) + ((deficit + 0xFFFFFFFFULL) >> 32);
        }
        const bool negative = k > un[j + n];
        un[j + n] = static_cast<u32>(un[j + n] - k);

        q[j] = static_cast<u32>(qhat);

        if (negative) {
            q[j] -= 1;
            u64 carry = 0;
            for (usize i = 0; i < n; ++i) {
                const u64 t = static_cast<u64>(un[i + j]) + vn[i] + carry;
                un[i + j] = static_cast<u32>(t);
                carry = t >> 32;
            }
            un[j + n] = static_cast<u32>(un[j + n] + carry);
        }
    }

    while (q.size() > 1 && q.back() == 0) q.pop_back();

    r.assign(n, u32{0});
    for (usize i = 0; i + 1 < n; ++i) {
        r[i] = s ? (un[i] >> s) | (un[i + 1] << (32 - s)) : un[i];
    }
    r[n - 1] = un[n - 1] >> s;
    while (r.size() > 1 && r.back() == 0) r.pop_back();
}

// Fixed xorshift stream for the above-2^64 probabilistic extension of is_prime.
u64 factor_rng_next(u64& state) noexcept {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 0x2545F4914F6CDD1DULL;
}

BigInt pow_mod_big(const BigInt& base, const BigInt& e, const BigInt& m) {
    BigInt r = BigInt(i64{1}) % m;
    BigInt b = base % m;
    const i32 bits = static_cast<i32>(e.bit_length());
    for (i32 i = bits - 1; i >= 0; --i) {
        r = r * r % m;
        if ((e >> static_cast<usize>(i)).low_limb() & 1ULL) {
            r = r * b % m;
        }
    }
    return r;
}

// Strong probable prime test of n (odd, > 2, a < n) to base a.
bool strong_probably_prime(const BigInt& n, const BigInt& a) {
    const BigInt one(i64{1});
    const BigInt n_minus_1 = n - one;
    BigInt d = n_minus_1;
    i32 s = 0;
    while ((d.low_limb() & 1ULL) == 0) {
        d = d >> 1;
        ++s;
    }
    BigInt x = pow_mod_big(a, d, n);
    if (x == one || x == n_minus_1) return true;
    for (i32 r = 1; r < s; ++r) {
        x = x * x % n;
        if (x == n_minus_1) return true;
    }
    return false;
}

// Pollard rho with Floyd cycles and batched gcd.
BigInt pollard_rho(const BigInt& n, const BigInt& c) {
    constexpr u64 kBatch = 96;
    const BigInt one(i64{1});
    BigInt x(i64{2});
    BigInt y(i64{2});
    for (;;) {
        const BigInt x0 = x;
        const BigInt y0 = y;
        BigInt q(i64{1});
        bool collided = false;
        for (u64 i = 0; i < kBatch; ++i) {
            x = (x * x + c) % n;
            y = (y * y + c) % n;
            y = (y * y + c) % n;
            const BigInt diff = ((x - y) % n).abs();
            if (diff.is_zero()) {
                collided = true;
                break;
            }
            q = q * diff % n;
        }
        if (!collided) {
            const BigInt g = BigInt::gcd(q, n);
            if (g == one) continue;
            if (g != n) return g;
        }
        x = x0;
        y = y0;
        for (u64 i = 0; i < kBatch; ++i) {
            x = (x * x + c) % n;
            y = (y * y + c) % n;
            y = (y * y + c) % n;
            const BigInt diff = ((x - y) % n).abs();
            if (diff.is_zero()) return n;
            const BigInt g = BigInt::gcd(diff, n);
            if (g != one) return g;
        }
    }
}

} // namespace

BigInt::BigInt(i64 v) noexcept {
    if (v == 0) return;
    m_negative = v < 0;
    m_limbs.push_back(m_negative ? (~static_cast<u64>(v) + 1ULL) : static_cast<u64>(v));
}

BigInt::BigInt(u64 v) noexcept {
    if (v != 0) m_limbs.push_back(v);
}

void BigInt::normalize() noexcept {
    while (!m_limbs.empty() && m_limbs.back() == 0) m_limbs.pop_back();
    if (m_limbs.empty()) m_negative = false;
}

i32 BigInt::signum() const noexcept {
    if (m_limbs.empty()) return 0;
    return m_negative ? -1 : 1;
}

usize BigInt::bit_length() const noexcept {
    if (m_limbs.empty()) return 0;
    return (m_limbs.size() - 1) * 64 + static_cast<usize>(64 - nlz64(m_limbs.back()));
}

i32 BigInt::compare_magnitude(const BigInt& o) const noexcept {
    if (m_limbs.size() != o.m_limbs.size()) {
        return m_limbs.size() < o.m_limbs.size() ? -1 : 1;
    }
    for (usize i = m_limbs.size(); i-- > 0;) {
        if (m_limbs[i] != o.m_limbs[i]) {
            return m_limbs[i] < o.m_limbs[i] ? -1 : 1;
        }
    }
    return 0;
}

BigInt BigInt::add_magnitude(const BigInt& a, const BigInt& b) {
    const LimbVec& big = a.m_limbs.size() >= b.m_limbs.size() ? a.m_limbs : b.m_limbs;
    const LimbVec& small = a.m_limbs.size() >= b.m_limbs.size() ? b.m_limbs : a.m_limbs;

    BigInt out;
    out.m_limbs.reserve(big.size() + 1);
    u64 carry = 0;
    for (usize i = 0; i < small.size(); ++i) {
        const u64 t = big[i] + small[i];
        const u64 c1 = (t < big[i]) ? 1ULL : 0ULL;
        const u64 t2 = t + carry;
        const u64 c2 = (t2 < t) ? 1ULL : 0ULL;
        out.m_limbs.push_back(t2);
        carry = c1 + c2;
    }
    for (usize i = small.size(); i < big.size(); ++i) {
        const u64 t = big[i] + carry;
        carry = (t < carry) ? 1ULL : 0ULL;
        out.m_limbs.push_back(t);
    }
    if (carry) out.m_limbs.push_back(1ULL);
    return out;
}

BigInt BigInt::sub_magnitude(const BigInt& a, const BigInt& b) {
    BigInt out;
    out.m_limbs.reserve(a.m_limbs.size());
    u64 borrow = 0;
    for (usize i = 0; i < b.m_limbs.size(); ++i) {
        const u64 t = a.m_limbs[i] - b.m_limbs[i];
        const u64 b1 = (a.m_limbs[i] < b.m_limbs[i]) ? 1ULL : 0ULL;
        const u64 t2 = t - borrow;
        const u64 b2 = (t < borrow) ? 1ULL : 0ULL;
        out.m_limbs.push_back(t2);
        borrow = b1 + b2;
    }
    for (usize i = b.m_limbs.size(); i < a.m_limbs.size(); ++i) {
        const u64 t = a.m_limbs[i] - borrow;
        borrow = (t > a.m_limbs[i]) ? 1ULL : 0ULL;
        out.m_limbs.push_back(t);
    }
    out.normalize();
    return out;
}

BigInt operator+(const BigInt& a, const BigInt& b) {
    BigInt out;
    if (a.m_negative == b.m_negative) {
        out = BigInt::add_magnitude(a, b);
        out.m_negative = a.m_negative;
        out.normalize();
        return out;
    }
    const i32 cmp = a.compare_magnitude(b);
    if (cmp == 0) return BigInt();
    if (cmp > 0) {
        out = BigInt::sub_magnitude(a, b);
        out.m_negative = a.m_negative;
    } else {
        out = BigInt::sub_magnitude(b, a);
        out.m_negative = b.m_negative;
    }
    out.normalize();
    return out;
}

BigInt operator-(const BigInt& a, const BigInt& b) {
    BigInt neg_b = b;
    neg_b.m_negative = !b.m_negative;
    neg_b.normalize();
    return a + neg_b;
}

BigInt operator*(const BigInt& a, const BigInt& b) {
    if (a.is_zero() || b.is_zero()) return BigInt();
    BigInt out;
    out.m_limbs.assign(a.m_limbs.size() + b.m_limbs.size(), 0);
    for (usize i = 0; i < a.m_limbs.size(); ++i) {
        u64 carry = 0;
        for (usize j = 0; j < b.m_limbs.size(); ++j) {
            u64 lo = 0;
            u64 hi = 0;
            mulu64(a.m_limbs[i], b.m_limbs[j], &lo, &hi);
            const u64 t = out.m_limbs[i + j] + lo;
            const u64 c1 = (t < lo) ? 1ULL : 0ULL;
            const u64 t2 = t + carry;
            const u64 c2 = (t2 < t) ? 1ULL : 0ULL;
            out.m_limbs[i + j] = t2;
            carry = hi + c1 + c2;
        }
        usize k = i + b.m_limbs.size();
        while (carry != 0) {
            const u64 t = out.m_limbs[k] + carry;
            carry = (t < carry) ? 1ULL : 0ULL;
            out.m_limbs[k] = t;
            ++k;
        }
    }
    out.m_negative = a.m_negative != b.m_negative;
    out.normalize();
    return out;
}

void BigInt::div_mod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r) {
    KEEL_ASSERT(!b.is_zero());
    q = BigInt();
    r = BigInt();
    if (a.is_zero()) return;

    if (a.m_limbs.size() == 1 && b.m_limbs.size() == 1) {
        const u64 ua = a.m_limbs[0];
        const u64 ub = b.m_limbs[0];
        const u64 uq = ua / ub;
        const u64 ur = ua % ub;
        if (uq != 0) {
            q.m_limbs.push_back(uq);
            q.m_negative = (a.m_negative != b.m_negative);
        }
        if (ur != 0) {
            r.m_limbs.push_back(ur);
            r.m_negative = a.m_negative;
        }
        return;
    }

    Array<u32> ud;
    Array<u32> vd;
    Array<u32> qd;
    Array<u32> rd;
    to_digits32(a.m_limbs, ud);
    to_digits32(b.m_limbs, vd);
    div_mod_digits(ud, vd, qd, rd);

    from_digits32(qd, q.m_limbs);
    from_digits32(rd, r.m_limbs);
    q.m_negative = (a.m_negative != b.m_negative) && !q.m_limbs.empty();
    r.m_negative = a.m_negative && !r.m_limbs.empty();
}

BigInt operator/(const BigInt& a, const BigInt& b) {
    BigInt q;
    BigInt r;
    BigInt::div_mod(a, b, q, r);
    return q;
}

BigInt operator%(const BigInt& a, const BigInt& b) {
    BigInt q;
    BigInt r;
    BigInt::div_mod(a, b, q, r);
    return r;
}

bool operator==(const BigInt& a, const BigInt& b) noexcept {
    return a.m_negative == b.m_negative && a.m_limbs == b.m_limbs;
}

bool operator!=(const BigInt& a, const BigInt& b) noexcept { return !(a == b); }

bool operator<(const BigInt& a, const BigInt& b) noexcept {
    if (a.m_negative != b.m_negative) return a.m_negative;
    const i32 cmp = a.compare_magnitude(b);
    return a.m_negative ? cmp > 0 : cmp < 0;
}

bool operator<=(const BigInt& a, const BigInt& b) noexcept { return !(b < a); }

bool operator>(const BigInt& a, const BigInt& b) noexcept { return b < a; }

bool operator>=(const BigInt& a, const BigInt& b) noexcept { return !(a < b); }

BigInt BigInt::abs() const {
    BigInt out = *this;
    out.m_negative = false;
    return out;
}

BigInt BigInt::operator-() const {
    BigInt out = *this;
    if (!out.m_limbs.empty()) out.m_negative = !out.m_negative;
    return out;
}

BigInt BigInt::operator<<(usize bits) const {
    if (is_zero()) return BigInt();
    const usize limb_off = bits / 64;
    const u32 bit_off = static_cast<u32>(bits % 64);
    BigInt out;
    out.m_negative = m_negative;
    out.m_limbs.assign(m_limbs.size() + limb_off + (bit_off ? 1 : 0), 0);
    if (bit_off == 0) {
        for (usize i = 0; i < m_limbs.size(); ++i) out.m_limbs[i + limb_off] = m_limbs[i];
    } else {
        u64 carry = 0;
        for (usize i = 0; i < m_limbs.size(); ++i) {
            out.m_limbs[i + limb_off] = (m_limbs[i] << bit_off) | carry;
            carry = m_limbs[i] >> (64 - bit_off);
        }
        out.m_limbs[m_limbs.size() + limb_off] = carry;
    }
    out.normalize();
    return out;
}

BigInt BigInt::operator>>(usize bits) const {
    if (is_zero()) return BigInt();
    const usize limb_off = bits / 64;
    const u32 bit_off = static_cast<u32>(bits % 64);

    bool dropped = false;
    for (usize i = 0; i < limb_off && i < m_limbs.size(); ++i) {
        if (m_limbs[i] != 0) dropped = true;
    }
    if (bit_off != 0 && limb_off < m_limbs.size()) {
        if ((m_limbs[limb_off] & ((1ULL << bit_off) - 1ULL)) != 0) dropped = true;
    }

    BigInt out;
    if (limb_off >= m_limbs.size()) {
        if (m_negative && dropped) return BigInt(i64{-1});
        return BigInt();
    }
    out.m_limbs.assign(m_limbs.size() - limb_off, 0);
    if (bit_off == 0) {
        for (usize i = 0; i < out.m_limbs.size(); ++i) out.m_limbs[i] = m_limbs[i + limb_off];
    } else {
        for (usize i = 0; i < out.m_limbs.size(); ++i) {
            const u64 hi =
                (i + limb_off + 1 < m_limbs.size()) ? (m_limbs[i + limb_off + 1] << (64 - bit_off)) : 0ULL;
            out.m_limbs[i] = (m_limbs[i + limb_off] >> bit_off) | hi;
        }
    }
    out.m_negative = m_negative;
    out.normalize();
    if (m_negative && dropped) {
        out = out - BigInt(i64{1});
    }
    return out;
}

BigInt BigInt::gcd(const BigInt& a, const BigInt& b) {
    if (a.m_limbs.size() <= 1 && b.m_limbs.size() <= 1) {
        u64 u = a.low_limb();
        u64 v = b.low_limb();
        while (v != 0) {
            const u64 t = u % v;
            u = v;
            v = t;
        }
        return BigInt(u);
    }
    BigInt x = a.abs();
    BigInt y = b.abs();
    while (!y.is_zero()) {
        if (x.m_limbs.size() <= 1 && y.m_limbs.size() <= 1) {
            u64 u = x.low_limb();
            u64 v = y.low_limb();
            while (v != 0) {
                const u64 t = u % v;
                u = v;
                v = t;
            }
            return BigInt(u);
        }
        BigInt q;
        BigInt r;
        div_mod(x, y, q, r);
        x = keel::move(y);
        y = keel::move(r);
    }
    return x;
}

BigInt BigInt::pow_u64(const BigInt& base, u64 exp) {
    BigInt result(i64{1});
    BigInt b = base;
    while (exp != 0) {
        if (exp & 1ULL) result = result * b;
        exp >>= 1;
        if (exp != 0) b = b * b;
    }
    return result;
}

BigInt BigInt::lcm(const BigInt& a, const BigInt& b) {
    if (a.is_zero() || b.is_zero()) return BigInt();
    return (a * b).abs() / gcd(a, b);
}

bool BigInt::sqrt_exact(BigInt& out) const {
    if (m_negative) return false;
    if (is_zero()) {
        out = BigInt();
        return true;
    }
    // Newton iteration on x -> (x + n/x)/2, started above sqrt(n).
    BigInt x = BigInt(i64{1}) << (bit_length() / 2 + 1);
    for (;;) {
        BigInt q;
        BigInt r;
        div_mod(*this, x, q, r);
        BigInt y = (x + q) >> 1;
        if (y >= x) break;
        x = keel::move(y);
    }
    if (x * x == *this) {
        out = keel::move(x);
        return true;
    }
    return false;
}

u64 BigInt::hash_bits() const noexcept {
    u64 h = 0x9E3779B97F4A7C15ULL;
    if (m_negative) h ^= 0x123456789ABCDEFULL;
    for (u64 l : m_limbs) h = (h ^ l) * 0x100000001B3ULL;
    return h;
}

bool BigInt::fits_i64(i64& out) const noexcept {
    if (m_limbs.empty()) {
        out = 0;
        return true;
    }
    if (m_limbs.size() > 1) return false;
    const u64 m = m_limbs[0];
    if (m_negative) {
        if (m > (1ULL << 63)) return false;
        out = (m == (1ULL << 63)) ? Limits<i64>::kMin : -static_cast<i64>(m);
        return true;
    }
    if (m > static_cast<u64>(Limits<i64>::kMax)) return false;
    out = static_cast<i64>(m);
    return true;
}

bool BigInt::fits_u64(u64& out) const noexcept {
    if (m_limbs.empty()) {
        out = 0;
        return true;
    }
    if (m_negative || m_limbs.size() > 1) return false;
    out = m_limbs[0];
    return true;
}

f64 BigInt::to_f64() const noexcept {
    if (m_limbs.empty()) return 0.0;
    const usize n = m_limbs.size();
    if (n == 1) {
        const f64 v = static_cast<f64>(m_limbs[0]);
        return m_negative ? -v : v;
    }
    const usize bits = bit_length();
    if (bits > 1024) {
        return m_negative ? -Limits<f64>::infinity() : Limits<f64>::infinity();
    }
    const u64 hi = m_limbs[n - 1];
    const u64 lo = m_limbs[n - 2];
    const i32 lz = nlz64(hi);
    // Align the top 64 bits of the magnitude into `top64`.
    const u64 top64 = (lz == 0) ? hi : ((hi << lz) | (lo >> (64 - lz)));
    // Round the top 64-bit word to 53 IEEE-754 mantissa bits (round-to-nearest-even).
    u64 mant53 = top64 >> 11;
    const u64 rem11 = top64 & 0x7FFULL;
    bool sticky = (rem11 & 0x3FFULL) != 0;
    if (!sticky && lz > 0 && (lo & ((1ULL << (64 - lz)) - 1ULL)) != 0) sticky = true;
    for (usize i = 0; !sticky && i + 2 < n; ++i) {
        if (m_limbs[i] != 0) sticky = true;
    }
    if ((rem11 & 0x400ULL) != 0 && (sticky || (mant53 & 1ULL) != 0)) ++mant53;
    const f64 scaled = keel::ldexp(static_cast<f64>(mant53), static_cast<i32>(bits) - 53);
    return m_negative ? -scaled : scaled;
}

String BigInt::to_string(Allocator* alloc) const {
    if (is_zero()) return String::from("0", alloc);
    const BigInt chunk(kPow10_19);
    Array<u64> chunks(alloc);
    BigInt m = abs();
    while (!m.is_zero()) {
        BigInt q;
        BigInt r;
        div_mod(m, chunk, q, r);
        chunks.push_back(r.low_limb());
        m = keel::move(q);
    }
    String digits(alloc);
    digits.reserve(chunks.size() * 19 + (m_negative ? 1 : 0));
    if (m_negative) digits.push_back('-');
    for (usize i = chunks.size(); i-- > 0;) {
        char buf[32];
        const usize n = keel::format_u64(buf, sizeof(buf), chunks[i]);
        if (i + 1 < chunks.size() && n < 19) {
            for (usize z = 0; z < 19 - n; ++z) digits.push_back('0');
        }
        digits.append(StringRef(buf, n));
    }
    return digits;
}

bool BigInt::try_parse(BigInt& out, StringRef s) {
    out = BigInt();
    usize i = 0;
    bool negative = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
        negative = (s[i] == '-');
        ++i;
    }
    if (i >= s.size()) return false;
    for (char ch : s.substr(i)) {
        if (ch < '0' || ch > '9') return false;
    }
    const usize ndigits = s.size() - i;
    usize lead = ndigits % 19;
    if (lead == 0) lead = 19;

    BigInt v;
    usize pos = i;
    while (pos < s.size()) {
        const usize len = (pos == i) ? lead : 19;
        u64 c = 0;
        for (usize k = 0; k < len; ++k) {
            c = c * 10ULL + static_cast<u64>(s[pos + k] - '0');
        }
        v = v * BigInt(kPow10[len]) + BigInt(c);
        pos += len;
    }
    v.m_negative = negative;
    v.normalize();
    out = keel::move(v);
    return true;
}

bool BigInt::is_prime(const BigInt& n) {
    static const u64 kBases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < BigInt(i64{2})) return false;
    for (u64 b : kBases) {
        const BigInt bb(b);
        if (n == bb) return true;
        BigInt q;
        BigInt r;
        div_mod(n, bb, q, r);
        if (r.is_zero()) return false;
    }
    u64 fits = 0;
    const bool below_2e64 = n.fits_u64(fits);
    for (u64 b : kBases) {
        if (!strong_probably_prime(n, BigInt(b))) return false;
    }
    if (below_2e64) return true;
    // Probabilistic extension: 38 more fixed-stream bases.
    u64 rng = 0x9E3779B97F4A7C15ULL;
    const BigInt span = n - BigInt(i64{2});
    for (i32 i = 0; i < 38; ++i) {
        const BigInt a = BigInt(factor_rng_next(rng)) % span + BigInt(i64{2});
        if (!strong_probably_prime(n, a)) return false;
    }
    return true;
}

Array<PrimeFactor> BigInt::factorize(const BigInt& n_in, Allocator* alloc) {
    KEEL_ASSERT(!(n_in < BigInt(i64{1})));
    Array<PrimeFactor> out(alloc);
    if (n_in == BigInt(i64{1})) return out;
    BigInt n = n_in;
    static const u64 kSmall[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : kSmall) {
        i32 e = 0;
        for (;;) {
            BigInt q;
            BigInt r;
            div_mod(n, BigInt(p), q, r);
            if (!r.is_zero()) break;
            n = keel::move(q);
            ++e;
        }
        if (e != 0) out.push_back(PrimeFactor{BigInt(p), e});
    }
    Array<BigInt> rest(alloc);
    if (n != BigInt(i64{1})) rest.push_back(keel::move(n));
    while (!rest.empty()) {
        const BigInt m = rest.back();
        rest.pop_back();
        if (m == BigInt(i64{1})) continue;
        if (is_prime(m)) {
            bool merged = false;
            for (PrimeFactor& pe : out) {
                if (pe.prime == m) {
                    ++pe.exponent;
                    merged = true;
                    break;
                }
            }
            if (!merged) out.push_back(PrimeFactor{m, 1});
            continue;
        }
        BigInt d(i64{1});
        for (u64 c = 1;; ++c) {
            d = pollard_rho(m, BigInt(c));
            if (d != BigInt(i64{1}) && d != m) break;
        }
        rest.push_back(d);
        rest.push_back(m / d);
    }
    keel::sort(out, [](const PrimeFactor& a, const PrimeFactor& b) { return a.prime < b.prime; });
    return out;
}

} // namespace keel::bignum
