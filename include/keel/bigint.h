#pragma once

#include "keel/array.h"
#include "keel/core.h"
#include "keel/memory.h"
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
// ============================================================================
class BigInt {
public:
    BigInt() noexcept = default;
    BigInt(i64 v) noexcept;
    BigInt(u64 v) noexcept; // non-negative

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsSigned<T> && !IsSame<RemoveCv<T>, i64>)
    BigInt(T v) noexcept : BigInt(static_cast<i64>(v)) {}

    template <typename T>
        requires(IsIntegral<T> && !IsSame<RemoveCv<T>, bool> && IsUnsigned<T> && !IsSame<RemoveCv<T>, u64>)
    BigInt(T v) noexcept : BigInt(static_cast<u64>(v)) {}

    // Strict decimal parse: optional sign, one or more digits, no whitespace.
    [[nodiscard]] static bool try_parse(BigInt& out, StringRef s);

    [[nodiscard]] bool is_zero() const noexcept { return m_limbs.empty(); }
    [[nodiscard]] bool is_inline() const noexcept { return m_limbs.is_inline(); }
    [[nodiscard]] i32 signum() const noexcept;
    [[nodiscard]] bool is_negative() const noexcept { return m_negative; }
    [[nodiscard]] usize bit_length() const noexcept;

    // Low 64-bit limb of the magnitude (returns 0 when zero).
    [[nodiscard]] u64 low_limb() const noexcept { return m_limbs.empty() ? 0 : m_limbs[0]; }

    [[nodiscard]] bool fits_i64(i64& out) const noexcept;
    [[nodiscard]] bool fits_u64(u64& out) const noexcept;

    // Nearest IEEE-754 double (+/-INFINITY on overflow past 2^1024).
    [[nodiscard]] f64 to_f64() const noexcept;

    [[nodiscard]] String to_string(Allocator* alloc = nullptr) const;

    [[nodiscard]] BigInt abs() const;
    [[nodiscard]] BigInt operator-() const;

    friend BigInt operator+(const BigInt& a, const BigInt& b);
    friend BigInt operator-(const BigInt& a, const BigInt& b);
    friend BigInt operator*(const BigInt& a, const BigInt& b);
    friend BigInt operator/(const BigInt& a, const BigInt& b);
    friend BigInt operator%(const BigInt& a, const BigInt& b);

    BigInt& operator+=(const BigInt& o) { return *this = *this + o; }
    BigInt& operator-=(const BigInt& o) { return *this = *this - o; }
    BigInt& operator*=(const BigInt& o) { return *this = *this * o; }
    BigInt& operator/=(const BigInt& o) { return *this = *this / o; }
    BigInt& operator%=(const BigInt& o) { return *this = *this % o; }

    friend bool operator==(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator!=(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator<(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator<=(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator>(const BigInt& a, const BigInt& b) noexcept;
    friend bool operator>=(const BigInt& a, const BigInt& b) noexcept;

    [[nodiscard]] BigInt operator<<(usize bits) const;
    [[nodiscard]] BigInt operator>>(usize bits) const; // floor semantics for negatives

    static void div_mod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r);
    [[nodiscard]] static BigInt gcd(const BigInt& a, const BigInt& b);         // non-negative result
    [[nodiscard]] static BigInt pow_u64(const BigInt& base, u64 exp);          // 0^0 == 1
    [[nodiscard]] static BigInt lcm(const BigInt& a, const BigInt& b);         // non-negative result

    // Miller-Rabin primality test. Deterministic for n < 2^64 (bases {2..37});
    // above 2^64 adds 38 fixed-stream bases (error < 4^-50).
    [[nodiscard]] static bool is_prime(const BigInt& n);

    // Canonical factorization of n >= 1 into ascending (prime, exponent) pairs.
    // Deterministic (Pollard-Brent rho with fixed c-sequence).
    [[nodiscard]] static Array<PrimeFactor> factorize(const BigInt& n, Allocator* alloc = nullptr);

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
                default_allocator()->free_array(m_heap, static_cast<usize>(m_cap));
            }
        }

        LimbVec(const LimbVec& o) : m_size(o.m_size) {
            if (o.m_size <= kInlineCap) {
                m_cap = static_cast<u32>(kInlineCap);
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.data()[i];
            } else {
                m_cap = o.m_size;
                m_heap = default_allocator()->alloc_array<u64>(static_cast<usize>(m_cap));
                for (u32 i = 0; i < o.m_size; ++i) m_heap[i] = o.m_heap[i];
            }
        }

        LimbVec(LimbVec&& o) noexcept : m_size(o.m_size), m_cap(o.m_cap) {
            if (o.is_inline()) {
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.m_inline[i];
            } else {
                m_heap = o.m_heap;
                o.m_heap = nullptr;
                o.m_cap = static_cast<u32>(kInlineCap);
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
            u64* fresh = default_allocator()->alloc_array<u64>(static_cast<usize>(o.m_size));
            const u64* src = o.data();
            for (u32 i = 0; i < o.m_size; ++i) fresh[i] = src[i];
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, static_cast<usize>(m_cap));
            }
            m_heap = fresh;
            m_size = o.m_size;
            m_cap = o.m_size;
            return *this;
        }

        LimbVec& operator=(LimbVec&& o) noexcept {
            if (this == &o) return *this;
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, static_cast<usize>(m_cap));
            }
            m_size = o.m_size;
            m_cap = o.m_cap;
            m_heap = nullptr;
            if (o.is_inline()) {
                for (u32 i = 0; i < o.m_size; ++i) m_inline[i] = o.m_inline[i];
            } else {
                m_heap = o.m_heap;
                o.m_heap = nullptr;
                o.m_cap = static_cast<u32>(kInlineCap);
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
            const usize doubled = static_cast<usize>(m_cap) * 2;
            const u32 target = static_cast<u32>(new_cap > doubled ? new_cap : doubled);
            u64* fresh = default_allocator()->alloc_array<u64>(static_cast<usize>(target));
            const u64* src = data();
            for (u32 i = 0; i < m_size; ++i) fresh[i] = src[i];
            if (!is_inline()) {
                default_allocator()->free_array(m_heap, static_cast<usize>(m_cap));
            }
            m_heap = fresh;
            m_cap = target;
        }

        void push_back(u64 v) {
            if (m_size == m_cap) reserve(static_cast<usize>(m_cap) * 2);
            data()[m_size++] = v;
        }

        void assign(usize count, u64 v) {
            if (count > m_cap) reserve(count);
            u64* dst = data();
            for (usize i = 0; i < count; ++i) dst[i] = v;
            m_size = static_cast<u32>(count);
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
        u32 m_cap = static_cast<u32>(kInlineCap);
    };

    void normalize() noexcept;
    [[nodiscard]] i32 compare_magnitude(const BigInt& o) const noexcept;
    [[nodiscard]] static BigInt add_magnitude(const BigInt& a, const BigInt& b);
    [[nodiscard]] static BigInt sub_magnitude(const BigInt& a, const BigInt& b); // |a| >= |b|

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
