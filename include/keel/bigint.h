#pragma once

#include "keel/array.h"
#include "keel/core.h"
#include "keel/memory.h"
#include "keel/optional.h"
#include "keel/ref.h"
#include "keel/string.h"

namespace keel::bignum {

struct PrimeFactor;

// ============================================================================
// BigInt — signed arbitrary-precision integer.
//
// Representation: signed magnitude, little-endian base-2^64 limbs with 2-limb
// (128-bit) Small Buffer Optimization (SBO) inline on the stack.
// Invariants: no trailing zero limbs; zero is the empty limb buffer with
// m_negative == false (no negative zero).
//
// Division/modulo follow C truncated semantics:
//   a == (a / b) * b + (a % b), |a % b| < |b|, sign(a % b) == sign(a) or 0.
// Division by zero is a caller contract violation (KEEL_ASSERT in debug).
//
// Bitwise operators (~, &, |, ^, >>) follow infinite two's-complement
// semantics matching Python int (e.g., ~x == -(x + 1), -5 >> 1 == -3).
// ============================================================================
class BigInt {
public:
    BigInt() noexcept = default;
    BigInt(i64 v) noexcept;
    BigInt(u64 v) noexcept; // non-negative

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsSigned<T> && !IsSame<RemoveCv<T>, i64>)
    BigInt(T v) noexcept : BigInt(as<i64>(v)) {}

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsUnsigned<T> && !IsSame<RemoveCv<T>, u64>)
    BigInt(T v) noexcept : BigInt(as<u64>(v)) {}

    // Strict decimal parse: optional sign, one or more digits, no whitespace.
    [[nodiscard]] static bool try_parse(BigInt& out, StringRef s);

    // Strict base-N parse for base in [2, 36]: optional sign, one or more digits
    // (0-9, a-z/A-Z), no prefix, no whitespace.
    [[nodiscard]] static bool try_parse_radix(BigInt& out, StringRef s, i32 base);

    // Strict hexadecimal parse: optional sign, optional "0x"/"0X" prefix, one or more hex digits.
    [[nodiscard]] static bool try_parse_hex(BigInt& out, StringRef s);

    // Unsigned binary byte deserialization (returns non-negative BigInt).
    [[nodiscard]] static BigInt from_bytes_be(ArrayRef<u8> bytes);
    [[nodiscard]] static BigInt from_bytes_le(ArrayRef<u8> bytes);

    [[nodiscard]] bool is_zero() const noexcept { return m_limbs.empty(); }
    [[nodiscard]] bool is_inline() const noexcept { return m_limbs.is_inline(); }
    [[nodiscard]] i32 signum() const noexcept;
    [[nodiscard]] bool is_negative() const noexcept { return m_negative; }
    [[nodiscard]] usize bit_length() const noexcept;
    [[nodiscard]] usize byte_length() const noexcept { return (bit_length() + 7) / 8; }

    // Low 64-bit limb of the magnitude (returns 0 when zero).
    [[nodiscard]] u64 low_limb() const noexcept { return m_limbs.empty() ? 0 : m_limbs[0]; }

    // Two's-complement bit inspection and mutation.
    [[nodiscard]] bool test_bit(usize bit_index) const noexcept;
    void set_bit(usize bit_index);
    void clear_bit(usize bit_index);
    void flip_bit(usize bit_index);

    // Count of trailing zero bits in |this| (returns 0 when zero).
    [[nodiscard]] usize trailing_zeros() const noexcept;

    // Count of set bits in the magnitude |this| (matches Python int.bit_count()).
    [[nodiscard]] usize popcount() const noexcept;

    [[nodiscard]] bool fits_i64(i64& out) const noexcept;
    [[nodiscard]] bool fits_u64(u64& out) const noexcept;

    // Nearest IEEE-754 double (+/-INFINITY on overflow past 2^1024).
    [[nodiscard]] f64 to_f64() const noexcept;

    [[nodiscard]] String to_string(Allocator* alloc = nullptr) const;
    [[nodiscard]] String to_string_radix(i32 base, bool uppercase = false, Allocator* alloc = nullptr) const;
    [[nodiscard]] String to_hex(bool uppercase = false, bool prefix = false, Allocator* alloc = nullptr) const;

    // Encodes the unsigned magnitude |this| in big-endian / little-endian byte order,
    // zero-padded to at least min_len bytes (returns empty array for 0 when min_len == 0).
    [[nodiscard]] Array<u8> to_bytes_be(usize min_len = 0, Allocator* alloc = nullptr) const;
    [[nodiscard]] Array<u8> to_bytes_le(usize min_len = 0, Allocator* alloc = nullptr) const;

    [[nodiscard]] BigInt abs() const;
    [[nodiscard]] BigInt operator-() const;
    [[nodiscard]] BigInt operator~() const; // two's-complement NOT: -(x + 1)

    friend BigInt operator+(const BigInt& a, const BigInt& b);
    friend BigInt operator-(const BigInt& a, const BigInt& b);
    friend BigInt operator*(const BigInt& a, const BigInt& b);
    friend BigInt operator/(const BigInt& a, const BigInt& b);
    friend BigInt operator%(const BigInt& a, const BigInt& b);

    friend BigInt operator&(const BigInt& a, const BigInt& b);
    friend BigInt operator|(const BigInt& a, const BigInt& b);
    friend BigInt operator^(const BigInt& a, const BigInt& b);

    BigInt& operator+=(const BigInt& o) { return *this = *this + o; }
    BigInt& operator-=(const BigInt& o) { return *this = *this - o; }
    BigInt& operator*=(const BigInt& o) { return *this = *this * o; }
    BigInt& operator/=(const BigInt& o) { return *this = *this / o; }
    BigInt& operator%=(const BigInt& o) { return *this = *this % o; }

    BigInt& operator&=(const BigInt& o) { return *this = *this & o; }
    BigInt& operator|=(const BigInt& o) { return *this = *this | o; }
    BigInt& operator^=(const BigInt& o) { return *this = *this ^ o; }
    BigInt& operator<<=(usize bits) { return *this = *this << bits; }
    BigInt& operator>>=(usize bits) { return *this = *this >> bits; }

    friend bool operator==(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator!=(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator<(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator<=(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator>(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator>=(const BigInt& a, const BigInt& b) noexcept;

    [[nodiscard]] BigInt operator<<(usize bits) const;
    [[nodiscard]] BigInt operator>>(usize bits) const; // floor semantics for negatives

    // Alias-safe truncated division with remainder.
    static void div_mod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r);

    // Canonical Euclidean modulo in [0, |m| - 1] (m != 0).
    [[nodiscard]] static BigInt mod_euclid(const BigInt& a, const BigInt& m);

    [[nodiscard]] static BigInt gcd(const BigInt& a, const BigInt& b);         // non-negative result
    [[nodiscard]] static BigInt lcm(const BigInt& a, const BigInt& b);         // non-negative result

    // Extended Euclidean algorithm: returns g = gcd(|a|, |b|) >= 0 with a*x + b*y == g.
    [[nodiscard]] static BigInt gcd_ext(const BigInt& a, const BigInt& b, BigInt& x, BigInt& y);

    // Modular multiplicative inverse in [0, m - 1]. Returns kNullOpt if m <= 0 or gcd(a, m) != 1.
    [[nodiscard]] static Optional<BigInt> inv_mod(const BigInt& a, const BigInt& m);

    // Modular exponentiation (base^exp mod mod) in [0, mod - 1].
    // Supports negative exp when base is invertible modulo mod. Returns kNullOpt if mod <= 0 or not invertible.
    [[nodiscard]] static Optional<BigInt> pow_mod(const BigInt& base, const BigInt& exp, const BigInt& mod);

    // Jacobi symbol (a / n) in {-1, 0, +1} for odd n > 0 (KEEL_ASSERT in debug).
    [[nodiscard]] static i32 jacobi(const BigInt& a, const BigInt& n);

    [[nodiscard]] static BigInt pow_u64(const BigInt& base, u64 exp);          // 0^0 == 1

    // Miller-Rabin primality test. Deterministic for n < 2^64 (bases {2..37});
    // above 2^64 adds 38 fixed-stream bases (error < 4^-50).
    [[nodiscard]] static bool is_prime(const BigInt& n);

    // Canonical factorization of n >= 1 into ascending (prime, exponent) pairs.
    // Deterministic (Pollard-Brent rho with fixed c-sequence).
    [[nodiscard]] static Array<PrimeFactor> factorize(const BigInt& n, Allocator* alloc = nullptr);

    // CSPRNG generators backed by keel::os_random_bytes.
    // random_bits(bits): uniform integer in [0, 2^bits - 1].
    // random_below(upper): unbiased uniform integer in [0, upper - 1] (kNullOpt if upper <= 0).
    // random_prime(bits): random probable prime with exact bit_length() == bits (kNullOpt if bits < 2).
    [[nodiscard]] static Optional<BigInt> random_bits(usize bits);
    [[nodiscard]] static Optional<BigInt> random_below(const BigInt& upper);
    [[nodiscard]] static Optional<BigInt> random_prime(usize bits);

    // Floor integer square root floor(sqrt(*this)). Returns kNullOpt for negatives.
    [[nodiscard]] Optional<BigInt> isqrt() const;

    // Exact integer square root (returns false for negatives and non-squares).
    [[nodiscard]] bool sqrt_exact(BigInt& out) const;

    // Content hash for hash tables / hash-consing.
    [[nodiscard]] u64 hash_bits() const noexcept;

private:
    // 2-limb (128-bit) Small Buffer Optimization container for base-2^64 limbs.
    class LimbVec {
    public:
        static constexpr usize kInlineCap = 2;

        LimbVec() noexcept = default;

        ~LimbVec() noexcept {
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, as<usize>(m_cap));
            }
        }

        LimbVec(const LimbVec& o) : m_size(o.m_size) {
            if (o.m_size <= kInlineCap) {
                m_cap = as<u32>(kInlineCap);
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.data()[i];
            } else {
                m_cap = o.m_size;
                m_heap = default_allocator()->alloc_array<u64>(as<usize>(m_cap));
                for (u32 i = 0; i < o.m_size; ++i) m_heap[i] = o.m_heap[i];
            }
        }

        LimbVec(LimbVec&& o) noexcept : m_size(o.m_size), m_cap(o.m_cap) {
            if (o.is_inline()) {
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.m_inline[i];
            } else {
                m_heap = o.m_heap;
                o.m_heap = nullptr;
                o.m_cap = as<u32>(kInlineCap);
            }
            o.m_size = 0;
        }

        LimbVec& operator=(const LimbVec& o) {
            if (this == &o) return *this;
            if (o.m_size <= m_cap) {
                u64* dst = data();
                const u64* src = o.data();
                for (u32 i = 0; i < o.m_size; ++i) dst[i] = src[i];
                m_size = o.m_size;
                return *this;
            }
            u64* fresh = default_allocator()->alloc_array<u64>(as<usize>(o.m_size));
            const u64* src = o.data();
            for (u32 i = 0; i < o.m_size; ++i) fresh[i] = src[i];
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, as<usize>(m_cap));
            }
            m_heap = fresh;
            m_size = o.m_size;
            m_cap = o.m_size;
            return *this;
        }

        LimbVec& operator=(LimbVec&& o) noexcept {
            if (this == &o) return *this;
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, as<usize>(m_cap));
            }
            m_size = o.m_size;
            m_cap = o.m_cap;
            m_heap = nullptr;
            if (o.is_inline()) {
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.m_inline[i];
            } else {
                m_heap = o.m_heap;
                o.m_heap = nullptr;
                o.m_cap = as<u32>(kInlineCap);
            }
            o.m_size = 0;
            return *this;
        }

        [[nodiscard]] bool is_inline() const noexcept { return m_cap <= kInlineCap; }
        [[nodiscard]] bool empty() const noexcept { return m_size == 0; }
        [[nodiscard]] usize size() const noexcept { return m_size; }

        [[nodiscard]] u64* data() noexcept { return is_inline() ? m_inline : m_heap; }
        [[nodiscard]] const u64* data() const noexcept { return is_inline() ? m_inline : m_heap; }

        [[nodiscard]] u64& operator[](usize i) noexcept {
            KEEL_BOUNDS_ASSERT(i < m_size);
            return data()[i];
        }
        [[nodiscard]] const u64& operator[](usize i) const noexcept {
            KEEL_BOUNDS_ASSERT(i < m_size);
            return data()[i];
        }

        [[nodiscard]] u64& back() noexcept {
            KEEL_BOUNDS_ASSERT(m_size > 0);
            return data()[m_size - 1];
        }
        [[nodiscard]] const u64& back() const noexcept {
            KEEL_BOUNDS_ASSERT(m_size > 0);
            return data()[m_size - 1];
        }

        [[nodiscard]] const u64* begin() const noexcept { return data(); }
        [[nodiscard]] const u64* end() const noexcept { return data() + m_size; }

        void clear() noexcept { m_size = 0; }

        void pop_back() noexcept {
            KEEL_BOUNDS_ASSERT(m_size > 0);
            --m_size;
        }

        void reserve(usize new_cap) {
            if (new_cap <= m_cap) return;
            KEEL_ASSERT(new_cap <= Limits<u32>::kMax);
            const usize doubled = as<usize>(m_cap) * 2;
            const usize chosen = new_cap > doubled ? new_cap : doubled;
            const u32 target = as<u32>(chosen > Limits<u32>::kMax ? as<usize>(Limits<u32>::kMax) : chosen);
            u64* fresh = default_allocator()->alloc_array<u64>(as<usize>(target));
            const u64* src = data();
            for (u32 i = 0; i < m_size; ++i) fresh[i] = src[i];
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, as<usize>(m_cap));
            }
            m_heap = fresh;
            m_cap = target;
        }

        void push_back(u64 v) {
            if (m_size == m_cap) reserve(as<usize>(m_cap) * 2);
            data()[m_size++] = v;
        }

        void assign(usize count, u64 v) {
            if (count > m_cap) reserve(count);
            u64* dst = data();
            for (usize i = 0; i < count; ++i) dst[i] = v;
            m_size = as<u32>(count);
        }

        void resize(usize count, u64 v = 0) {
            if (count > m_cap) reserve(count);
            u64* dst = data();
            for (usize i = m_size; i < count; ++i) dst[i] = v;
            m_size = as<u32>(count);
        }

        [[nodiscard]] bool operator==(const LimbVec& o) const noexcept {
            if (m_size != o.m_size) return false;
            const u64* a = data();
            const u64* b = o.data();
            for (u32 i = 0; i < m_size; ++i) {
                if (a[i] != b[i]) return false;
            }
            return true;
        }

    private:
        u64 m_inline[kInlineCap] = {0, 0};
        u64* m_heap = nullptr;
        u32 m_size = 0;
        u32 m_cap = as<u32>(kInlineCap);
    };

    void normalize() noexcept;
    [[nodiscard]] i32 compare_magnitude(const BigInt& o) const noexcept;
    [[nodiscard]] static BigInt add_magnitude(const BigInt& a, const BigInt& b);
    [[nodiscard]] static BigInt sub_magnitude(const BigInt& a, const BigInt& b); // |a| >= |b|
    [[nodiscard]] static BigInt mul_magnitude(const BigInt& a, const BigInt& b);
    [[nodiscard]] static BigInt slice_limbs(const BigInt& a, usize start, usize count);
    [[nodiscard]] static u32 div_mod_u32_mag(const BigInt& a, u32 d, BigInt& q);
    [[nodiscard]] static u32 mod_u32_mag(const BigInt& a, u32 d) noexcept;
    [[nodiscard]] static BigInt bitwise_op(const BigInt& a, const BigInt& b, i32 kind);

    LimbVec m_limbs;
    bool m_negative = false;
};

struct PrimeFactor {
    BigInt prime;
    i32 exponent = 0;

    [[nodiscard]] bool operator==(const PrimeFactor& o) const noexcept {
        return prime == o.prime && exponent == o.exponent;
    }
    [[nodiscard]] bool operator!=(const PrimeFactor& o) const noexcept { return !(*this == o); }
};

} // namespace keel::bignum
