#include <keel/algo.h>
#include <keel/bigint.h>
#include <keel/bits.h>
#include <keel/core.h>
#include <keel/math.h>
#include <keel/number.h>
#include <keel/os_random.h>

namespace keel::bignum {

namespace {

constexpr usize kKaratsubaThreshold = 32;

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

constexpr char kDigitsLower[36] = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
    'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
};

constexpr char kDigitsUpper[36] = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
    'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
};

i32 decode_digit(char ch) noexcept {
    if (ch >= '0' && ch <= '9') return as<i32>(ch - '0');
    if (ch >= 'a' && ch <= 'z') return 10 + as<i32>(ch - 'a');
    if (ch >= 'A' && ch <= 'Z') return 10 + as<i32>(ch - 'A');
    return -1;
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
        d.push_back(as<u32>(l));
        d.push_back(as<u32>(l >> 32));
    }
    while (!d.empty() && d.back() == 0) d.pop_back();
}

template <typename LimbContainer>
void from_digits32(const Array<u32>& d, LimbContainer& limbs) {
    limbs.assign((d.size() + 1) / 2, 0);
    for (usize i = 0; i < d.size(); ++i) {
        if (i & 1) {
            limbs[i >> 1] |= as<u64>(d[i]) << 32;
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
            q[j] = as<u32>(t / v[0]);
            k = t % v[0];
        }
        r.assign(1, as<u32>(k));
        while (q.size() > 1 && q.back() == 0) q.pop_back();
        return;
    }

    const i32 s = keel::countl_zero(v[n - 1]);
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
        const u64 num = (as<u64>(un[j + n]) << 32) | un[j + n - 1];
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
            un[i + j] = as<u32>(cur - sub);
            const u64 deficit = sub > cur ? sub - cur : 0;
            k = (p >> 32) + ((deficit + 0xFFFFFFFFULL) >> 32);
        }
        const bool negative = k > un[j + n];
        un[j + n] = as<u32>(un[j + n] - k);

        q[j] = as<u32>(qhat);

        if (negative) {
            q[j] -= 1;
            u64 carry = 0;
            for (usize i = 0; i < n; ++i) {
                const u64 t = as<u64>(un[i + j]) + vn[i] + carry;
                un[i + j] = as<u32>(t);
                carry = t >> 32;
            }
            un[j + n] = as<u32>(un[j + n] + carry);
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

// Fast modular exponentiation for non-negative exponent e and modulus m > 1.
BigInt pow_mod_pos(const BigInt& base, const BigInt& e, const BigInt& m) {
    BigInt b = BigInt::mod_euclid(base, m);
    if (e.is_zero()) return BigInt(i64{1});
    if (b.is_zero()) return BigInt(i64{0});

    const usize bits = e.bit_length();
    if (bits < 32) {
        BigInt r(i64{1});
        for (usize i = bits; i-- > 0;) {
            r = (r * r) % m;
            if (e.test_bit(i)) {
                r = (r * b) % m;
            }
        }
        return r;
    }

    // 4-bit fixed window exponentiation.
    BigInt table[16];
    table[0] = BigInt(i64{1});
    table[1] = b;
    for (usize k = 2; k < 16; ++k) {
        table[k] = (table[k - 1] * b) % m;
    }

    BigInt r(i64{1});
    const usize top_nib = (bits - 1) / 4;
    for (usize nib = top_nib + 1; nib-- > 0;) {
        if (nib != top_nib) {
            r = (r * r) % m;
            r = (r * r) % m;
            r = (r * r) % m;
            r = (r * r) % m;
        }
        const usize base_bit = nib * 4;
        u32 w = 0;
        if (e.test_bit(base_bit + 3)) w |= 8u;
        if (e.test_bit(base_bit + 2)) w |= 4u;
        if (e.test_bit(base_bit + 1)) w |= 2u;
        if (e.test_bit(base_bit)) w |= 1u;
        if (w != 0u) {
            r = (r * table[w]) % m;
        }
    }
    return r;
}

// Strong probable prime test of n (odd, > 2, a < n) to base a.
bool strong_probably_prime(const BigInt& n, const BigInt& a) {
    const BigInt one(i64{1});
    const BigInt n_minus_1 = n - one;
    const usize s = n_minus_1.trailing_zeros();
    const BigInt d = n_minus_1 >> s;
    BigInt x = pow_mod_pos(a, d, n);
    if (x == one || x == n_minus_1) return true;
    for (usize r = 1; r < s; ++r) {
        x = (x * x) % n;
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
    m_limbs.push_back(m_negative ? (~as<u64>(v) + 1ULL) : as<u64>(v));
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
    return (m_limbs.size() - 1) * 64 + as<usize>(64 - keel::countl_zero(m_limbs.back()));
}

usize BigInt::trailing_zeros() const noexcept {
    if (m_limbs.empty()) return 0;
    for (usize i = 0; i < m_limbs.size(); ++i) {
        if (m_limbs[i] != 0) {
            return i * 64 + as<usize>(keel::countr_zero(m_limbs[i]));
        }
    }
    return 0;
}

usize BigInt::popcount() const noexcept {
    usize total = 0;
    for (u64 l : m_limbs) {
        total += as<usize>(keel::popcount(l));
    }
    return total;
}

bool BigInt::test_bit(usize bit_index) const noexcept {
    const usize limb = bit_index / 64;
    const u32 off = as<u32>(bit_index % 64);
    if (!m_negative) {
        if (limb >= m_limbs.size()) return false;
        return ((m_limbs[limb] >> off) & 1ULL) != 0;
    }
    const usize tz = trailing_zeros();
    if (bit_index < tz) return false;
    if (bit_index == tz) return true;
    const bool mag_bit = (limb < m_limbs.size()) && (((m_limbs[limb] >> off) & 1ULL) != 0);
    return !mag_bit;
}

void BigInt::set_bit(usize bit_index) {
    if (!m_negative) {
        const usize limb = bit_index / 64;
        const u32 off = as<u32>(bit_index % 64);
        if (limb >= m_limbs.size()) m_limbs.resize(limb + 1, 0);
        m_limbs[limb] |= (1ULL << off);
        return;
    }
    *this |= (BigInt(u64{1}) << bit_index);
}

void BigInt::clear_bit(usize bit_index) {
    if (!m_negative) {
        const usize limb = bit_index / 64;
        const u32 off = as<u32>(bit_index % 64);
        if (limb >= m_limbs.size()) return;
        m_limbs[limb] &= ~(1ULL << off);
        normalize();
        return;
    }
    *this &= ~(BigInt(u64{1}) << bit_index);
}

void BigInt::flip_bit(usize bit_index) {
    if (!m_negative) {
        const usize limb = bit_index / 64;
        const u32 off = as<u32>(bit_index % 64);
        if (limb >= m_limbs.size()) m_limbs.resize(limb + 1, 0);
        m_limbs[limb] ^= (1ULL << off);
        normalize();
        return;
    }
    *this ^= (BigInt(u64{1}) << bit_index);
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

BigInt BigInt::slice_limbs(const BigInt& a, usize start, usize count) {
    BigInt out;
    if (start >= a.m_limbs.size() || count == 0) return out;
    const usize rem = a.m_limbs.size() - start;
    const usize n = count < rem ? count : rem;
    out.m_limbs.reserve(n);
    for (usize i = 0; i < n; ++i) {
        out.m_limbs.push_back(a.m_limbs[start + i]);
    }
    out.normalize();
    return out;
}

BigInt BigInt::mul_magnitude(const BigInt& a, const BigInt& b) {
    const usize na = a.m_limbs.size();
    const usize nb = b.m_limbs.size();
    if (na == 0 || nb == 0) return BigInt();

    if (na < kKaratsubaThreshold || nb < kKaratsubaThreshold) {
        BigInt out;
        out.m_limbs.assign(na + nb, 0);
        for (usize i = 0; i < na; ++i) {
            u64 carry = 0;
            for (usize j = 0; j < nb; ++j) {
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
            usize k = i + nb;
            while (carry != 0) {
                const u64 t = out.m_limbs[k] + carry;
                carry = (t < carry) ? 1ULL : 0ULL;
                out.m_limbs[k] = t;
                ++k;
            }
        }
        out.normalize();
        return out;
    }

    const usize max_n = na > nb ? na : nb;
    const usize m = max_n / 2;

    // Unbalanced split: avoid redundant Karatsuba work when one operand fits in m limbs.
    if (m >= nb) {
        const BigInt a0 = slice_limbs(a, 0, m);
        const BigInt a1 = slice_limbs(a, m, na - m);
        const BigInt lo = mul_magnitude(a0, b);
        const BigInt hi = mul_magnitude(a1, b);
        return add_magnitude(lo, hi << (m * 64));
    }
    if (m >= na) {
        const BigInt b0 = slice_limbs(b, 0, m);
        const BigInt b1 = slice_limbs(b, m, nb - m);
        const BigInt lo = mul_magnitude(a, b0);
        const BigInt hi = mul_magnitude(a, b1);
        return add_magnitude(lo, hi << (m * 64));
    }

    const BigInt a0 = slice_limbs(a, 0, m);
    const BigInt a1 = slice_limbs(a, m, na - m);
    const BigInt b0 = slice_limbs(b, 0, m);
    const BigInt b1 = slice_limbs(b, m, nb - m);

    const BigInt z0 = mul_magnitude(a0, b0);
    const BigInt z2 = mul_magnitude(a1, b1);
    const BigInt s_a = add_magnitude(a0, a1);
    const BigInt s_b = add_magnitude(b0, b1);
    const BigInt z1 = mul_magnitude(s_a, s_b);
    const BigInt mid = sub_magnitude(sub_magnitude(z1, z2), z0);

    BigInt res = add_magnitude(z0, mid << (m * 64));
    res = add_magnitude(res, z2 << (m * 128));
    res.normalize();
    return res;
}

u32 BigInt::div_mod_u32_mag(const BigInt& a, u32 d, BigInt& q) {
    KEEL_ASSERT(d != 0u);
    BigInt q_out;
    if (a.m_limbs.empty()) {
        q = BigInt();
        return 0u;
    }
    q_out.m_limbs.assign(a.m_limbs.size(), 0);
    u64 rem = 0;
    const u64 div = as<u64>(d);
    for (usize i = a.m_limbs.size(); i-- > 0;) {
        const u64 limb = a.m_limbs[i];
        const u64 hi_cur = (rem << 32) | (limb >> 32);
        const u64 q_hi = hi_cur / div;
        rem = hi_cur % div;
        const u64 lo_cur = (rem << 32) | (limb & 0xFFFFFFFFULL);
        const u64 q_lo = lo_cur / div;
        rem = lo_cur % div;
        q_out.m_limbs[i] = (q_hi << 32) | q_lo;
    }
    q_out.normalize();
    q = keel::move(q_out);
    return as<u32>(rem);
}

u32 BigInt::mod_u32_mag(const BigInt& a, u32 d) noexcept {
    KEEL_ASSERT(d != 0u);
    u64 rem = 0;
    const u64 div = as<u64>(d);
    for (usize i = a.m_limbs.size(); i-- > 0;) {
        const u64 limb = a.m_limbs[i];
        rem = ((rem << 32) | (limb >> 32)) % div;
        rem = ((rem << 32) | (limb & 0xFFFFFFFFULL)) % div;
    }
    return as<u32>(rem);
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
    if (b.is_zero()) return a;
    if (a.is_zero()) return -b;
    const bool neg_b = !b.m_negative;
    BigInt out;
    if (a.m_negative == neg_b) {
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
        out.m_negative = neg_b;
    }
    out.normalize();
    return out;
}

BigInt operator*(const BigInt& a, const BigInt& b) {
    if (a.is_zero() || b.is_zero()) return BigInt();
    BigInt out = BigInt::mul_magnitude(a, b);
    out.m_negative = a.m_negative != b.m_negative;
    out.normalize();
    return out;
}

void BigInt::div_mod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r) {
    KEEL_ASSERT(!b.is_zero());
    if (a.is_zero()) {
        q = BigInt();
        r = BigInt();
        return;
    }

    BigInt q_out;
    BigInt r_out;
    const bool q_neg = (a.m_negative != b.m_negative);
    const bool r_neg = a.m_negative;

    if (a.m_limbs.size() == 1 && b.m_limbs.size() == 1) {
        const u64 ua = a.m_limbs[0];
        const u64 ub = b.m_limbs[0];
        const u64 uq = ua / ub;
        const u64 ur = ua % ub;
        if (uq != 0) {
            q_out.m_limbs.push_back(uq);
            q_out.m_negative = q_neg;
        }
        if (ur != 0) {
            r_out.m_limbs.push_back(ur);
            r_out.m_negative = r_neg;
        }
        q = keel::move(q_out);
        r = keel::move(r_out);
        return;
    }

    Array<u32> ud;
    Array<u32> vd;
    Array<u32> qd;
    Array<u32> rd;
    to_digits32(a.m_limbs, ud);
    to_digits32(b.m_limbs, vd);
    div_mod_digits(ud, vd, qd, rd);

    from_digits32(qd, q_out.m_limbs);
    from_digits32(rd, r_out.m_limbs);
    q_out.m_negative = q_neg && !q_out.m_limbs.empty();
    r_out.m_negative = r_neg && !r_out.m_limbs.empty();
    q = keel::move(q_out);
    r = keel::move(r_out);
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

BigInt BigInt::mod_euclid(const BigInt& a, const BigInt& m) {
    KEEL_ASSERT(!m.is_zero());
    const BigInt abs_m = m.abs();
    BigInt r = a % abs_m;
    if (r.is_negative()) {
        r += abs_m;
    }
    return r;
}

BigInt BigInt::operator~() const {
    // Two's-complement NOT: ~x == -(x + 1).
    return -(*this + BigInt(i64{1}));
}

// kind: 0 = AND, 1 = OR, 2 = XOR (exact infinite two's-complement semantics).
BigInt BigInt::bitwise_op(const BigInt& a, const BigInt& b, i32 kind) {
    const usize na = a.m_limbs.size();
    const usize nb = b.m_limbs.size();

    // Fast path for non-negative operands.
    if (!a.m_negative && !b.m_negative) {
        BigInt out;
        if (kind == 0) {
            const usize n = na < nb ? na : nb;
            out.m_limbs.reserve(n);
            for (usize i = 0; i < n; ++i) {
                out.m_limbs.push_back(a.m_limbs[i] & b.m_limbs[i]);
            }
        } else if (kind == 1) {
            const usize n = na > nb ? na : nb;
            out.m_limbs.reserve(n);
            for (usize i = 0; i < n; ++i) {
                const u64 la = i < na ? a.m_limbs[i] : 0ULL;
                const u64 lb = i < nb ? b.m_limbs[i] : 0ULL;
                out.m_limbs.push_back(la | lb);
            }
        } else {
            const usize n = na > nb ? na : nb;
            out.m_limbs.reserve(n);
            for (usize i = 0; i < n; ++i) {
                const u64 la = i < na ? a.m_limbs[i] : 0ULL;
                const u64 lb = i < nb ? b.m_limbs[i] : 0ULL;
                out.m_limbs.push_back(la ^ lb);
            }
        }
        out.normalize();
        return out;
    }

    // General infinite two's-complement path.
    const usize n = na > nb ? na : nb;
    const u64 sign_a = a.m_negative ? 1ULL : 0ULL;
    const u64 sign_b = b.m_negative ? 1ULL : 0ULL;
    u64 sign_res = 0ULL;
    if (kind == 0) {
        sign_res = sign_a & sign_b;
    } else if (kind == 1) {
        sign_res = sign_a | sign_b;
    } else {
        sign_res = sign_a ^ sign_b;
    }

    u64 carry_a = a.m_negative ? 1ULL : 0ULL;
    u64 carry_b = b.m_negative ? 1ULL : 0ULL;
    u64 carry_out = sign_res ? 1ULL : 0ULL;

    BigInt out;
    out.m_limbs.reserve(n + 1);
    for (usize i = 0; i < n; ++i) {
        u64 ta = i < na ? a.m_limbs[i] : 0ULL;
        if (a.m_negative) {
            const u64 inv = ~ta;
            ta = inv + carry_a;
            carry_a = (ta < inv) ? 1ULL : 0ULL;
        }
        u64 tb = i < nb ? b.m_limbs[i] : 0ULL;
        if (b.m_negative) {
            const u64 inv = ~tb;
            tb = inv + carry_b;
            carry_b = (tb < inv) ? 1ULL : 0ULL;
        }
        u64 tr = 0ULL;
        if (kind == 0) {
            tr = ta & tb;
        } else if (kind == 1) {
            tr = ta | tb;
        } else {
            tr = ta ^ tb;
        }
        if (sign_res != 0ULL) {
            const u64 inv = ~tr;
            const u64 mag = inv + carry_out;
            carry_out = (mag < inv) ? 1ULL : 0ULL;
            out.m_limbs.push_back(mag);
        } else {
            out.m_limbs.push_back(tr);
        }
    }
    if (sign_res != 0ULL && carry_out != 0ULL) {
        out.m_limbs.push_back(1ULL);
    }
    out.m_negative = (sign_res != 0ULL);
    out.normalize();
    return out;
}

BigInt operator&(const BigInt& a, const BigInt& b) { return BigInt::bitwise_op(a, b, 0); }
BigInt operator|(const BigInt& a, const BigInt& b) { return BigInt::bitwise_op(a, b, 1); }
BigInt operator^(const BigInt& a, const BigInt& b) { return BigInt::bitwise_op(a, b, 2); }

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
    if (bits == 0) return *this;
    const usize limb_off = bits / 64;
    const u32 bit_off = as<u32>(bits % 64);
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
    if (bits == 0) return *this;
    const usize limb_off = bits / 64;
    const u32 bit_off = as<u32>(bits % 64);

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

BigInt BigInt::gcd_ext(const BigInt& a, const BigInt& b, BigInt& x, BigInt& y) {
    const bool neg_a = a.m_negative;
    const bool neg_b = b.m_negative;
    BigInt old_r = a.abs();
    BigInt r = b.abs();
    BigInt old_s(i64{1});
    BigInt s(i64{0});
    BigInt old_t(i64{0});
    BigInt t(i64{1});

    while (!r.is_zero()) {
        BigInt q;
        BigInt rem;
        div_mod(old_r, r, q, rem);
        old_r = keel::move(r);
        r = keel::move(rem);

        BigInt next_s = old_s - q * s;
        old_s = keel::move(s);
        s = keel::move(next_s);

        BigInt next_t = old_t - q * t;
        old_t = keel::move(t);
        t = keel::move(next_t);
    }

    if (neg_a) old_s = -old_s;
    if (neg_b) old_t = -old_t;
    x = keel::move(old_s);
    y = keel::move(old_t);
    return old_r;
}

Optional<BigInt> BigInt::inv_mod(const BigInt& a, const BigInt& m) {
    if (m.is_zero() || m.is_negative()) return kNullOpt;
    if (m == BigInt(i64{1})) return BigInt(i64{0});
    BigInt x;
    BigInt y;
    const BigInt g = gcd_ext(mod_euclid(a, m), m, x, y);
    if (g != BigInt(i64{1})) return kNullOpt;
    return mod_euclid(x, m);
}

Optional<BigInt> BigInt::pow_mod(const BigInt& base, const BigInt& exp, const BigInt& mod) {
    if (mod.is_zero() || mod.is_negative()) return kNullOpt;
    if (mod == BigInt(i64{1})) return BigInt(i64{0});
    if (exp.is_negative()) {
        Optional<BigInt> inv = inv_mod(base, mod);
        if (!inv) return kNullOpt;
        return pow_mod_pos(*inv, -exp, mod);
    }
    return pow_mod_pos(base, exp, mod);
}

i32 BigInt::jacobi(const BigInt& a, const BigInt& n) {
    KEEL_ASSERT(!n.is_negative() && !n.is_zero() && (n.low_limb() & 1ULL) == 1ULL);
    if (n.is_negative() || n.is_zero() || (n.low_limb() & 1ULL) == 0ULL) return 0;
    if (n == BigInt(i64{1})) return 1;

    BigInt x = mod_euclid(a, n);
    BigInt m = n;
    i32 result = 1;

    while (!x.is_zero()) {
        const usize tz = x.trailing_zeros();
        if (tz > 0) {
            x >>= tz;
            if ((tz & 1u) != 0u) {
                const u64 m_mod8 = m.low_limb() & 7ULL;
                if (m_mod8 == 3ULL || m_mod8 == 5ULL) result = -result;
            }
        }
        BigInt tmp = keel::move(x);
        x = keel::move(m);
        m = keel::move(tmp);
        if ((x.low_limb() & 3ULL) == 3ULL && (m.low_limb() & 3ULL) == 3ULL) {
            result = -result;
        }
        x = x % m;
    }

    return (m == BigInt(i64{1})) ? result : 0;
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
    return (a.abs() / gcd(a, b)) * b.abs();
}

Optional<BigInt> BigInt::isqrt() const {
    if (m_negative) return kNullOpt;
    if (is_zero()) return BigInt();
    // Newton iteration on x -> (x + n/x)/2, started strictly above sqrt(n).
    BigInt x = BigInt(i64{1}) << (bit_length() / 2 + 1);
    for (;;) {
        BigInt q;
        BigInt r;
        div_mod(*this, x, q, r);
        BigInt y = (x + q) >> 1;
        if (y >= x) break;
        x = keel::move(y);
    }
    return x;
}

bool BigInt::sqrt_exact(BigInt& out) const {
    Optional<BigInt> root = isqrt();
    if (!root) return false;
    if ((*root) * (*root) == *this) {
        out = keel::move(*root);
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
        out = (m == (1ULL << 63)) ? Limits<i64>::kMin : -as<i64>(m);
        return true;
    }
    if (m > as<u64>(Limits<i64>::kMax)) return false;
    out = as<i64>(m);
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
        const f64 v = as<f64>(m_limbs[0]);
        return m_negative ? -v : v;
    }
    const usize bits = bit_length();
    if (bits > 1024) {
        return m_negative ? -Limits<f64>::infinity() : Limits<f64>::infinity();
    }
    const u64 hi = m_limbs[n - 1];
    const u64 lo = m_limbs[n - 2];
    const i32 lz = keel::countl_zero(hi);
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
    const f64 scaled = keel::ldexp(as<f64>(mant53), as<i32>(bits) - 53);
    return m_negative ? -scaled : scaled;
}

BigInt BigInt::from_bytes_be(ArrayRef<u8> bytes) {
    usize start = 0;
    while (start < bytes.size() && bytes[start] == 0u) ++start;
    const usize len = bytes.size() - start;
    if (len == 0) return BigInt();

    const usize nlimbs = (len + 7) / 8;
    BigInt out;
    out.m_limbs.assign(nlimbs, 0ULL);
    for (usize i = 0; i < len; ++i) {
        const u8 b = bytes[bytes.size() - 1 - i];
        out.m_limbs[i / 8] |= as<u64>(b) << ((i % 8) * 8);
    }
    out.normalize();
    return out;
}

BigInt BigInt::from_bytes_le(ArrayRef<u8> bytes) {
    usize len = bytes.size();
    while (len > 0 && bytes[len - 1] == 0u) --len;
    if (len == 0) return BigInt();

    const usize nlimbs = (len + 7) / 8;
    BigInt out;
    out.m_limbs.assign(nlimbs, 0ULL);
    for (usize i = 0; i < len; ++i) {
        out.m_limbs[i / 8] |= as<u64>(bytes[i]) << ((i % 8) * 8);
    }
    out.normalize();
    return out;
}

Array<u8> BigInt::to_bytes_be(usize min_len, Allocator* alloc) const {
    const usize blen = byte_length();
    const usize out_len = blen > min_len ? blen : min_len;
    Array<u8> out(out_len, u8{0}, alloc);
    for (usize i = 0; i < blen; ++i) {
        out[out_len - 1 - i] = as<u8>((m_limbs[i / 8] >> ((i % 8) * 8)) & 0xFFULL);
    }
    return out;
}

Array<u8> BigInt::to_bytes_le(usize min_len, Allocator* alloc) const {
    const usize blen = byte_length();
    const usize out_len = blen > min_len ? blen : min_len;
    Array<u8> out(out_len, u8{0}, alloc);
    for (usize i = 0; i < blen; ++i) {
        out[i] = as<u8>((m_limbs[i / 8] >> ((i % 8) * 8)) & 0xFFULL);
    }
    return out;
}

String BigInt::to_string(Allocator* alloc) const {
    if (is_zero()) return String::from("0", alloc);
    constexpr u32 kChunkBase = 1000000000u; // 10^9 < 2^32
    Array<u32> chunks(alloc);
    BigInt m = abs();
    while (!m.is_zero()) {
        BigInt q;
        const u32 rem = div_mod_u32_mag(m, kChunkBase, q);
        chunks.push_back(rem);
        m = keel::move(q);
    }
    String digits(alloc);
    digits.reserve(chunks.size() * 9 + (m_negative ? 1 : 0));
    if (m_negative) digits.push_back('-');
    for (usize i = chunks.size(); i-- > 0;) {
        char buf[16];
        const usize n = keel::format_u64(buf, sizeof(buf), as<u64>(chunks[i]));
        if (i + 1 < chunks.size() && n < 9) {
            for (usize z = 0; z < 9 - n; ++z) digits.push_back('0');
        }
        digits.append(StringRef(buf, n));
    }
    return digits;
}

String BigInt::to_string_radix(i32 base, bool uppercase, Allocator* alloc) const {
    KEEL_ASSERT(base >= 2 && base <= 36);
    if (base < 2 || base > 36) return String(alloc);
    if (base == 10) return to_string(alloc);
    if (is_zero()) return String::from("0", alloc);

    const char* alphabet = uppercase ? kDigitsUpper : kDigitsLower;
    const u32 ubase = as<u32>(base);

    // Direct O(N) limb extraction for power-of-2 bases (2, 4, 8, 16, 32).
    if ((ubase & (ubase - 1u)) == 0u) {
        const u32 bits_per_digit = as<u32>(keel::countr_zero(ubase));
        const u64 mask = (1ULL << bits_per_digit) - 1ULL;
        const usize total_bits = bit_length();
        const usize ndigits = (total_bits + bits_per_digit - 1) / bits_per_digit;
        String out(alloc);
        out.reserve(ndigits + (m_negative ? 1 : 0));
        if (m_negative) out.push_back('-');
        for (usize d = ndigits; d-- > 0;) {
            const usize bit_pos = d * bits_per_digit;
            const usize limb = bit_pos / 64;
            const u32 off = as<u32>(bit_pos % 64);
            u64 val = m_limbs[limb] >> off;
            if (off + bits_per_digit > 64 && limb + 1 < m_limbs.size()) {
                val |= m_limbs[limb + 1] << (64 - off);
            }
            out.push_back(alphabet[val & mask]);
        }
        return out;
    }

    // General base: divide repeatedly by largest base^k fitting in u32.
    u32 chunk_base = 1u;
    usize chunk_digits = 0;
    while (chunk_base <= Limits<u32>::kMax / ubase) {
        chunk_base *= ubase;
        ++chunk_digits;
    }

    Array<u32> chunks(alloc);
    BigInt m = abs();
    while (!m.is_zero()) {
        BigInt q;
        chunks.push_back(div_mod_u32_mag(m, chunk_base, q));
        m = keel::move(q);
    }

    String out(alloc);
    out.reserve(chunks.size() * chunk_digits + (m_negative ? 1 : 0));
    if (m_negative) out.push_back('-');

    for (usize i = chunks.size(); i-- > 0;) {
        char buf[32];
        u32 v = chunks[i];
        usize len = 0;
        do {
            buf[len++] = alphabet[v % ubase];
            v /= ubase;
        } while (v != 0u);
        if (i + 1 < chunks.size()) {
            while (len < chunk_digits) {
                buf[len++] = '0';
            }
        }
        for (usize k = len; k-- > 0;) {
            out.push_back(buf[k]);
        }
    }
    return out;
}

String BigInt::to_hex(bool uppercase, bool prefix, Allocator* alloc) const {
    if (!prefix) return to_string_radix(16, uppercase, alloc);
    String mag = abs().to_string_radix(16, uppercase, alloc);
    String out(alloc);
    out.reserve(mag.size() + (m_negative ? 3 : 2));
    if (m_negative) out.push_back('-');
    out.append("0x");
    out.append(mag);
    return out;
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
            c = c * 10ULL + as<u64>(s[pos + k] - '0');
        }
        v = v * BigInt(kPow10[len]) + BigInt(c);
        pos += len;
    }
    v.m_negative = negative;
    v.normalize();
    out = keel::move(v);
    return true;
}

bool BigInt::try_parse_radix(BigInt& out, StringRef s, i32 base) {
    out = BigInt();
    if (base < 2 || base > 36) return false;
    if (base == 10) return try_parse(out, s);

    usize i = 0;
    bool negative = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
        negative = (s[i] == '-');
        ++i;
    }
    if (i >= s.size()) return false;
    for (char ch : s.substr(i)) {
        const i32 d = decode_digit(ch);
        if (d < 0 || d >= base) return false;
    }

    const u64 ubase = as<u64>(base);
    u64 max_pow = 1ULL;
    usize chunk_max = 0;
    while (max_pow <= Limits<u64>::kMax / ubase) {
        max_pow *= ubase;
        ++chunk_max;
    }

    const usize ndigits = s.size() - i;
    usize lead = ndigits % chunk_max;
    if (lead == 0) lead = chunk_max;

    BigInt v;
    usize pos = i;
    while (pos < s.size()) {
        const usize len = (pos == i) ? lead : chunk_max;
        u64 c = 0;
        u64 p = 1ULL;
        for (usize k = 0; k < len; ++k) {
            c = c * ubase + as<u64>(decode_digit(s[pos + k]));
            p *= ubase;
        }
        v = v * BigInt(p) + BigInt(c);
        pos += len;
    }
    v.m_negative = negative;
    v.normalize();
    out = keel::move(v);
    return true;
}

bool BigInt::try_parse_hex(BigInt& out, StringRef s) {
    out = BigInt();
    usize i = 0;
    bool negative = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
        negative = (s[i] == '-');
        ++i;
    }
    if (i + 2 <= s.size() && s[i] == '0' && (s[i + 1] == 'x' || s[i + 1] == 'X')) {
        i += 2;
    }
    if (i >= s.size()) return false;
    // Reject embedded signs after 0x prefix.
    if (s[i] == '+' || s[i] == '-') return false;

    BigInt mag;
    if (!try_parse_radix(mag, s.substr(i), 16)) return false;
    mag.m_negative = negative;
    mag.normalize();
    out = keel::move(mag);
    return true;
}

bool BigInt::is_prime(const BigInt& n) {
    static const u32 kBases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < BigInt(i64{2})) return false;
    for (u32 b : kBases) {
        if (n.m_limbs.size() == 1 && n.m_limbs[0] == as<u64>(b)) return true;
        if (mod_u32_mag(n, b) == 0u) return false;
    }
    u64 fits = 0;
    const bool below_2e64 = n.fits_u64(fits);
    for (u32 b : kBases) {
        if (!strong_probably_prime(n, BigInt(as<u64>(b)))) return false;
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

Optional<BigInt> BigInt::random_bits(usize bits) {
    if (bits == 0) return BigInt();
    const usize nbytes = (bits + 7) / 8;
    Array<u8> buf(nbytes, u8{0});
    if (keel::os_random_bytes(MutArrayRef<u8>(buf)).is_err()) {
        return kNullOpt;
    }
    const usize extra = nbytes * 8 - bits;
    if (extra > 0) {
        buf[0] = as<u8>(buf[0] & as<u8>(0xFFu >> extra));
    }
    return from_bytes_be(ArrayRef<u8>(buf));
}

Optional<BigInt> BigInt::random_below(const BigInt& upper) {
    if (upper.is_zero() || upper.is_negative()) return kNullOpt;
    if (upper == BigInt(i64{1})) return BigInt();
    const usize bits = (upper - BigInt(i64{1})).bit_length();
    for (;;) {
        Optional<BigInt> cand = random_bits(bits);
        if (!cand) return kNullOpt;
        if (*cand < upper) return cand;
    }
}

Optional<BigInt> BigInt::random_prime(usize bits) {
    if (bits < 2) return kNullOpt;
    if (bits == 2) {
        u8 b = 0;
        if (keel::os_random_bytes(MutArrayRef<u8>(&b, 1)).is_err()) return kNullOpt;
        return BigInt(i64{(b & 1u) ? 3 : 2});
    }

    static const u32 kSmallPrimes[] = {
        3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,
        71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
        163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
    };

    for (;;) {
        Optional<BigInt> cand = random_bits(bits);
        if (!cand) return kNullOpt;
        cand->set_bit(bits - 1);
        cand->set_bit(0);

        bool composite = false;
        for (u32 p : kSmallPrimes) {
            if (cand->m_limbs.size() == 1 && cand->m_limbs[0] == as<u64>(p)) {
                return cand;
            }
            if (mod_u32_mag(*cand, p) == 0u) {
                composite = true;
                break;
            }
        }
        if (composite) continue;
        if (is_prime(*cand)) return cand;
    }
}

Array<PrimeFactor> BigInt::factorize(const BigInt& n_in, Allocator* alloc) {
    KEEL_ASSERT(!(n_in < BigInt(i64{1})));
    Array<PrimeFactor> out(alloc);
    if (n_in == BigInt(i64{1})) return out;
    BigInt n = n_in;
    static const u32 kSmall[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u32 p : kSmall) {
        i32 e = 0;
        for (;;) {
            BigInt q;
            const u32 rem = div_mod_u32_mag(n, p, q);
            if (rem != 0u) break;
            n = keel::move(q);
            ++e;
        }
        if (e != 0) out.push_back(PrimeFactor{BigInt(as<u64>(p)), e});
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
