#include <keel/bignum.h>
#include <keel/format.h>
#include <keel/fs.h>
#include <keel/math.h>
#include <keel/number.h>

namespace {

keel::i32 g_failed = 0;

void check(bool cond, const char* expr, keel::i32 line) {
    if (!cond) {
        keel::eprintln("FAIL [line {}]: {}", line, expr);
        ++g_failed;
    }
}

#define CHECK(cond) check((cond), #cond, __LINE__)

keel::u64 f64_to_bits(keel::f64 v) noexcept {
    keel::u64 bits = 0;
    keel::mem_copy(&bits, &v, sizeof(bits));
    return bits;
}

keel::f64 bits_to_f64(keel::u64 bits) noexcept {
    keel::f64 v = 0.0;
    keel::mem_copy(&v, &bits, sizeof(v));
    return v;
}

keel::String bytes_to_hex(keel::ArrayRef<keel::u8> data) {
    static constexpr char kHex[16] = {'0', '1', '2', '3', '4', '5', '6', '7',
                                      '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
    keel::String out;
    out.reserve(data.size() * 2);
    for (keel::u8 b : data) {
        out.push_back(kHex[(b >> 4) & 0xFu]);
        out.push_back(kHex[b & 0xFu]);
    }
    return out;
}

bool hex_to_bytes(keel::StringRef hex, keel::Array<keel::u8>& out) {
    out.clear();
    if (hex == "-") return true; // sentinel for empty byte slice in batch CLI
    if ((hex.size() & 1u) != 0u) return false;
    out.reserve(hex.size() / 2);
    auto nib = [](char c) noexcept -> keel::i32 {
        if (c >= '0' && c <= '9') return keel::as<keel::i32>(c - '0');
        if (c >= 'a' && c <= 'f') return 10 + keel::as<keel::i32>(c - 'a');
        if (c >= 'A' && c <= 'F') return 10 + keel::as<keel::i32>(c - 'A');
        return -1;
    };
    for (keel::usize i = 0; i < hex.size(); i += 2) {
        const keel::i32 hi = nib(hex[i]);
        const keel::i32 lo = nib(hex[i + 1]);
        if (hi < 0 || lo < 0) return false;
        out.push_back(keel::as<keel::u8>((hi << 4) | lo));
    }
    return true;
}

keel::String u64_to_hex16(keel::u64 v) {
    static constexpr char kHex[16] = {'0', '1', '2', '3', '4', '5', '6', '7',
                                      '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
    keel::String out;
    out.reserve(16);
    for (keel::i32 s = 60; s >= 0; s -= 4) {
        out.push_back(kHex[(v >> keel::as<keel::u32>(s)) & 0xFULL]);
    }
    return out;
}

bool parse_u64_hex(keel::StringRef s, keel::u64& out) {
    out = 0;
    if (s.empty() || s.size() > 16) return false;
    for (char c : s) {
        keel::u64 d = 0;
        if (c >= '0' && c <= '9') {
            d = keel::as<keel::u64>(c - '0');
        } else if (c >= 'a' && c <= 'f') {
            d = 10ULL + keel::as<keel::u64>(c - 'a');
        } else if (c >= 'A' && c <= 'F') {
            d = 10ULL + keel::as<keel::u64>(c - 'A');
        } else {
            return false;
        }
        out = (out << 4) | d;
    }
    return true;
}

void test_bigint_basic_and_sbo() {
    using keel::bignum::BigInt;

    BigInt zero;
    CHECK(zero.is_zero());
    CHECK(zero.is_inline());
    CHECK(zero.signum() == 0);
    CHECK(zero.bit_length() == 0);
    CHECK(zero.byte_length() == 0);
    CHECK(zero.trailing_zeros() == 0);
    CHECK(zero.popcount() == 0);
    CHECK(zero.to_string() == "0");

    BigInt one(1);
    BigInt neg_min(keel::Limits<keel::i64>::kMin);
    CHECK(one.is_inline());
    CHECK(neg_min.is_inline());
    CHECK(neg_min.is_negative());
    CHECK(neg_min.to_string() == "-9223372036854775808");
    CHECK(neg_min.trailing_zeros() == 63);
    CHECK(neg_min.popcount() == 1);

    keel::i64 i_out = 0;
    CHECK(neg_min.fits_i64(i_out) && i_out == keel::Limits<keel::i64>::kMin);
    keel::u64 u_out = 0;
    CHECK(!neg_min.fits_u64(u_out));

    // 2 limbs inline, 3 limbs heap.
    BigInt two_limbs = (BigInt(1) << 127) - BigInt(1);
    CHECK(two_limbs.is_inline());
    CHECK(two_limbs.bit_length() == 127);
    CHECK(two_limbs.byte_length() == 16);
    CHECK(two_limbs.popcount() == 127);

    BigInt three_limbs = BigInt(1) << 128;
    CHECK(!three_limbs.is_inline());
    CHECK(three_limbs.bit_length() == 129);
    CHECK(three_limbs.byte_length() == 17);
    CHECK(three_limbs.trailing_zeros() == 128);
    CHECK(three_limbs.popcount() == 1);

    BigInt shrunk = three_limbs >> 64;
    BigInt copy_shrunk = shrunk;
    CHECK(copy_shrunk.is_inline());
    CHECK(copy_shrunk == (BigInt(1) << 64));
}

void test_bigint_arith_div_shift_karatsuba() {
    using keel::bignum::BigInt;

    BigInt a;
    BigInt b;
    CHECK(BigInt::try_parse(a, "1234567890123456789012345678901234567890"));
    CHECK(BigInt::try_parse(b, "-987654321098765432109876543210"));
    CHECK(a.to_string() == "1234567890123456789012345678901234567890");
    CHECK(b.to_string() == "-987654321098765432109876543210");

    BigInt q;
    BigInt r;
    BigInt::div_mod(a, b, q, r);
    CHECK(q * b + r == a);
    CHECK(r.abs() < b.abs());
    CHECK(!r.is_negative());

    // Alias-safety of div_mod when output aliases input.
    BigInt aliased_a = a;
    BigInt aliased_b = b;
    BigInt::div_mod(aliased_a, aliased_b, aliased_a, aliased_b);
    CHECK(aliased_a == q && aliased_b == r);

    BigInt::div_mod(-a, b, q, r);
    CHECK(q * b + r == -a);
    CHECK(r.is_negative() || r.is_zero());

    // Euclidean modulo.
    CHECK(BigInt::mod_euclid(BigInt(-1), BigInt(5)) == BigInt(4));
    CHECK(BigInt::mod_euclid(BigInt(-5), BigInt(5)) == BigInt(0));
    CHECK(BigInt::mod_euclid(BigInt(7), BigInt(-5)) == BigInt(2));
    CHECK(BigInt::mod_euclid(BigInt(-7), BigInt(-5)) == BigInt(3));

    // Knuth D borrow edge case: high-bit limbs.
    BigInt u = (BigInt(1) << 192) - BigInt(1);
    BigInt v = (BigInt(1) << 96) + BigInt(1);
    BigInt::div_mod(u, v, q, r);
    CHECK(q * v + r == u);
    CHECK(q == (BigInt(1) << 96) - BigInt(1));
    CHECK(r.is_zero());

    // Karatsuba path (>= 32 limbs = 2048 bits): verify (2^2500 - 1) * (2^2500 + 1) == 2^5000 - 1.
    BigInt km = (BigInt(1) << 2500) - BigInt(1);
    BigInt kp = (BigInt(1) << 2500) + BigInt(1);
    BigInt kprod = km * kp;
    CHECK(kprod == (BigInt(1) << 5000) - BigInt(1));
    BigInt::div_mod(kprod, kp, q, r);
    CHECK(q == km && r.is_zero());

    // Unbalanced Karatsuba (4000-bit * 2100-bit): verify via division identity.
    BigInt u_big = ((BigInt(1) << 4000) - BigInt(12345)) * BigInt(-1);
    BigInt v_big = (BigInt(1) << 2100) + BigInt(67891);
    BigInt uv = u_big * v_big;
    BigInt::div_mod(uv, v_big, q, r);
    CHECK(q == u_big && r.is_zero());

    // Negative right shift uses floor semantics.
    CHECK((BigInt(-1) >> 1) == BigInt(-1));
    CHECK((BigInt(-5) >> 1) == BigInt(-3));
    CHECK((BigInt(-4) >> 1) == BigInt(-2));

    // GCD / LCM / pow / isqrt / sqrt_exact.
    CHECK(BigInt::gcd(BigInt(48), BigInt(-18)) == BigInt(6));
    CHECK(BigInt::lcm(BigInt(12), BigInt(18)) == BigInt(36));
    CHECK(BigInt::pow_u64(BigInt(0), 0) == BigInt(1));
    CHECK(BigInt::pow_u64(BigInt(3), 40).to_string() == "12157665459056928801");

    BigInt sq_root;
    BigInt sq = a * a;
    CHECK(sq.sqrt_exact(sq_root) && sq_root == a);
    CHECK(!(sq + BigInt(1)).sqrt_exact(sq_root));
    CHECK(!BigInt(-4).sqrt_exact(sq_root));

    auto isq1 = (sq + a).isqrt();
    CHECK(isq1.has_value() && *isq1 == a);
    auto isq2 = (sq - BigInt(1)).isqrt();
    CHECK(isq2.has_value() && *isq2 == a - BigInt(1));
    CHECK(!BigInt(-1).isqrt().has_value());
    CHECK(BigInt(0).isqrt().value() == BigInt(0));

    // Invalid parse inputs.
    BigInt bad;
    CHECK(!BigInt::try_parse(bad, ""));
    CHECK(!BigInt::try_parse(bad, "-"));
    CHECK(!BigInt::try_parse(bad, " 123"));
    CHECK(!BigInt::try_parse(bad, "12a"));
}

void test_bigint_bitwise() {
    using keel::bignum::BigInt;

    // Positive bitwise ops.
    BigInt x(0xF0F0ULL);
    BigInt y(0x0FF0ULL);
    CHECK((x & y) == BigInt(0x00F0ULL));
    CHECK((x | y) == BigInt(0xFFF0ULL));
    CHECK((x ^ y) == BigInt(0xFF00ULL));
    CHECK(~BigInt(0) == BigInt(-1));
    CHECK(~BigInt(-1) == BigInt(0));
    CHECK(~BigInt(42) == BigInt(-43));
    CHECK(~BigInt(-43) == BigInt(42));

    // Negative two's-complement bitwise ops (matching Python).
    CHECK((BigInt(-5) & BigInt(-3)) == BigInt(-7));
    CHECK((BigInt(-5) | BigInt(-3)) == BigInt(-1));
    CHECK((BigInt(-5) ^ BigInt(-3)) == BigInt(6));
    CHECK((BigInt(-5) & BigInt(3)) == BigInt(3));
    CHECK((BigInt(-5) | BigInt(3)) == BigInt(-5));
    CHECK((BigInt(-5) ^ BigInt(3)) == BigInt(-8));

    // Multi-limb borrow/carry boundary: -(2^64 - 1) & -2 == -(2^64).
    BigInt neg_max64 = -((BigInt(1) << 64) - BigInt(1));
    CHECK((neg_max64 & BigInt(-2)) == -(BigInt(1) << 64));

    // Bit inspection / mutation on positive and negative values.
    BigInt m(0);
    m.set_bit(0);
    m.set_bit(63);
    m.set_bit(128);
    CHECK(m.test_bit(0) && !m.test_bit(1) && m.test_bit(63) && m.test_bit(128) && !m.test_bit(129));
    CHECK(m.popcount() == 3);
    m.clear_bit(128);
    CHECK(!m.test_bit(128) && m.bit_length() == 64);
    m.flip_bit(63);
    CHECK(m == BigInt(1));

    BigInt neg5(-5); // ...11111011
    CHECK(neg5.test_bit(0));
    CHECK(neg5.test_bit(1));
    CHECK(!neg5.test_bit(2));
    CHECK(neg5.test_bit(3));
    CHECK(neg5.test_bit(200));
    neg5.set_bit(2);
    CHECK(neg5 == BigInt(-1));
    neg5.clear_bit(0);
    CHECK(neg5 == BigInt(-2));
    neg5.flip_bit(1);
    CHECK(neg5 == BigInt(-4));

    // Compound operators.
    BigInt c(15);
    c &= BigInt(6);
    CHECK(c == BigInt(6));
    c |= BigInt(8);
    CHECK(c == BigInt(14));
    c ^= BigInt(3);
    CHECK(c == BigInt(13));
    c <<= 4;
    CHECK(c == BigInt(208));
    c >>= 2;
    CHECK(c == BigInt(52));
}

void test_bigint_radix_hex_and_bytes() {
    using keel::bignum::BigInt;

    // Hex parsing & formatting.
    BigInt h;
    CHECK(BigInt::try_parse_hex(h, "0xdeadBEEF0123456789abcdef"));
    CHECK(h.to_hex(false, false) == "deadbeef0123456789abcdef");
    CHECK(h.to_hex(true, true) == "0xDEADBEEF0123456789ABCDEF");
    CHECK(BigInt::try_parse_hex(h, "-0X1a2B"));
    CHECK(h == BigInt(-0x1a2b));
    CHECK(h.to_hex(false, true) == "-0x1a2b");
    CHECK(BigInt(0).to_hex(false, true) == "0x0");
    CHECK(BigInt::try_parse_hex(h, "+ff") && h == BigInt(255));
    CHECK(BigInt::try_parse_hex(h, "-0x0") && h.is_zero() && !h.is_negative());
    CHECK(BigInt::try_parse_hex(h, "+0X0") && h.is_zero() && !h.is_negative());
    CHECK(!BigInt::try_parse_hex(h, "0x"));
    CHECK(!BigInt::try_parse_hex(h, "-0x"));
    CHECK(!BigInt::try_parse_hex(h, "0x-1"));
    CHECK(!BigInt::try_parse_hex(h, "0x1g"));

    // General radix 2..36.
    BigInt v;
    CHECK(BigInt::try_parse_radix(v, "1011010111001", 2));
    CHECK(v.to_string_radix(2) == "1011010111001");
    CHECK(v.to_string_radix(8) == "13271");
    CHECK(v.to_string_radix(16) == "16b9");
    CHECK(BigInt::try_parse_radix(v, "-KeEl36Z", 36));
    CHECK(v.to_string_radix(36, true) == "-KEEL36Z");
    CHECK(!BigInt::try_parse_radix(v, "12", 1));
    CHECK(!BigInt::try_parse_radix(v, "12", 37));
    CHECK(!BigInt::try_parse_radix(v, "12", 2));

    // Binary byte serialization (BE & LE + zero padding).
    const keel::u8 raw_be[] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09};
    BigInt from_be = BigInt::from_bytes_be(keel::ArrayRef<keel::u8>(raw_be));
    CHECK(from_be.to_hex() == "10203040506070809");
    CHECK(from_be.byte_length() == 9);
    auto out_be = from_be.to_bytes_be();
    CHECK(out_be.size() == 9 && out_be[0] == 0x01 && out_be[8] == 0x09);
    auto padded_be = from_be.to_bytes_be(16);
    CHECK(padded_be.size() == 16 && padded_be[0] == 0x00 && padded_be[7] == 0x01 && padded_be[15] == 0x09);

    auto out_le = from_be.to_bytes_le(12);
    CHECK(out_le.size() == 12 && out_le[0] == 0x09 && out_le[8] == 0x01 && out_le[11] == 0x00);
    CHECK(BigInt::from_bytes_le(keel::ArrayRef<keel::u8>(out_le)) == from_be);

    // Zero bytes.
    CHECK(BigInt::from_bytes_be(keel::ArrayRef<keel::u8>()).is_zero());
    CHECK(BigInt(0).to_bytes_be().empty());
    auto z4 = BigInt(0).to_bytes_be(4);
    CHECK(z4.size() == 4 && z4[0] == 0 && z4[3] == 0);
}

void test_bigint_modular_and_number_theory() {
    using keel::bignum::BigInt;

    // Extended GCD.
    BigInt x;
    BigInt y;
    BigInt g = BigInt::gcd_ext(BigInt(240), BigInt(46), x, y);
    CHECK(g == BigInt(2));
    CHECK(BigInt(240) * x + BigInt(46) * y == g);

    g = BigInt::gcd_ext(BigInt(-240), BigInt(46), x, y);
    CHECK(g == BigInt(2));
    CHECK(BigInt(-240) * x + BigInt(46) * y == g);

    // Modular inverse.
    auto inv1 = BigInt::inv_mod(BigInt(3), BigInt(11));
    CHECK(inv1.has_value() && *inv1 == BigInt(4));
    auto inv_neg = BigInt::inv_mod(BigInt(-3), BigInt(11));
    CHECK(inv_neg.has_value() && *inv_neg == BigInt(7));
    CHECK(!BigInt::inv_mod(BigInt(6), BigInt(9)).has_value());
    CHECK(!BigInt::inv_mod(BigInt(3), BigInt(0)).has_value());
    CHECK(!BigInt::inv_mod(BigInt(3), BigInt(-5)).has_value());
    CHECK(BigInt::inv_mod(BigInt(5), BigInt(1)).value() == BigInt(0));

    // Modular exponentiation (positive, zero, negative exponents, RSA-style roundtrip).
    CHECK(BigInt::pow_mod(BigInt(2), BigInt(10), BigInt(1000)).value() == BigInt(24));
    CHECK(BigInt::pow_mod(BigInt(3), BigInt(-1), BigInt(11)).value() == BigInt(4));
    CHECK(BigInt::pow_mod(BigInt(3), BigInt(-2), BigInt(11)).value() == BigInt(5));
    CHECK(!BigInt::pow_mod(BigInt(2), BigInt(-1), BigInt(4)).has_value());
    CHECK(!BigInt::pow_mod(BigInt(2), BigInt(3), BigInt(0)).has_value());
    CHECK(!BigInt::pow_mod(BigInt(2), BigInt(3), BigInt(-7)).has_value());
    CHECK(BigInt::pow_mod(BigInt(999), BigInt(0), BigInt(1)).value() == BigInt(0));

    // Jacobi symbol.
    CHECK(BigInt::jacobi(BigInt(1001), BigInt(9907)) == -1);
    CHECK(BigInt::jacobi(BigInt(19), BigInt(45)) == 1);
    CHECK(BigInt::jacobi(BigInt(6), BigInt(15)) == 0);
    CHECK(BigInt::jacobi(BigInt(0), BigInt(1)) == 1);
    CHECK(BigInt::jacobi(BigInt(-1), BigInt(7)) == -1);
    CHECK(BigInt::jacobi(BigInt(-1), BigInt(5)) == 1);
}

void test_bigint_primality_factorize_and_csprng() {
    using keel::bignum::BigInt;

    CHECK(!BigInt::is_prime(BigInt(-7)));
    CHECK(!BigInt::is_prime(BigInt(0)));
    CHECK(!BigInt::is_prime(BigInt(1)));
    CHECK(BigInt::is_prime(BigInt(2)));
    CHECK(BigInt::is_prime(BigInt(3)));
    CHECK(!BigInt::is_prime(BigInt(561)));                   // Carmichael number
    CHECK(!BigInt::is_prime(BigInt(341550071728321ULL)));    // strong pseudoprime to 2..17
    CHECK(BigInt::is_prime(BigInt(2305843009213693951ULL))); // M61 = 2^61 - 1

    CHECK(BigInt::factorize(BigInt(1)).empty());

    // 2^64 - 1 = 3 * 5 * 17 * 257 * 641 * 65537 * 6700417
    BigInt f64 = (BigInt(1) << 64) - BigInt(1);
    auto fac = BigInt::factorize(f64);
    CHECK(fac.size() == 7);
    BigInt prod(1);
    for (keel::usize i = 0; i < fac.size(); ++i) {
        CHECK(BigInt::is_prime(fac[i].prime));
        CHECK(fac[i].exponent == 1);
        if (i > 0) CHECK(fac[i - 1].prime < fac[i].prime);
        prod *= BigInt::pow_u64(fac[i].prime, keel::as<keel::u64>(fac[i].exponent));
    }
    CHECK(prod == f64);

    // 100-bit semiprime above 2^64.
    BigInt p1(1000000007ULL);
    BigInt p2(1000000009ULL);
    BigInt p3(1000000000039ULL);
    BigInt n = p1 * p1 * p2 * p3;
    auto fac2 = BigInt::factorize(n);
    CHECK(fac2.size() == 3);
    CHECK(fac2[0].prime == p1 && fac2[0].exponent == 2);
    CHECK(fac2[1].prime == p2 && fac2[1].exponent == 1);
    CHECK(fac2[2].prime == p3 && fac2[2].exponent == 1);

    // CSPRNG generators + mini RSA-256 keygen/encrypt/decrypt roundtrip!
    CHECK(BigInt::random_bits(0).value().is_zero());
    auto r128 = BigInt::random_bits(128);
    CHECK(r128.has_value() && r128->bit_length() <= 128 && !r128->is_negative());

    CHECK(!BigInt::random_below(BigInt(0)).has_value());
    CHECK(!BigInt::random_below(BigInt(-5)).has_value());
    CHECK(BigInt::random_below(BigInt(1)).value().is_zero());
    BigInt bound = (BigInt(1) << 100) + BigInt(12345);
    for (keel::i32 i = 0; i < 16; ++i) {
        auto rv = BigInt::random_below(bound);
        CHECK(rv.has_value() && !rv->is_negative() && *rv < bound);
    }

    CHECK(!BigInt::random_prime(0).has_value());
    CHECK(!BigInt::random_prime(1).has_value());
    auto rp2 = BigInt::random_prime(2);
    CHECK(rp2.has_value() && (*rp2 == BigInt(2) || *rp2 == BigInt(3)));

    auto rp_p = BigInt::random_prime(128);
    auto rp_q = BigInt::random_prime(128);
    CHECK(rp_p.has_value() && rp_p->bit_length() == 128 && BigInt::is_prime(*rp_p));
    CHECK(rp_q.has_value() && rp_q->bit_length() == 128 && BigInt::is_prime(*rp_q));
    if (rp_p && rp_q && *rp_p != *rp_q) {
        BigInt mod_n = (*rp_p) * (*rp_q);
        BigInt phi = (*rp_p - BigInt(1)) * (*rp_q - BigInt(1));
        BigInt pub_e(65537);
        auto priv_d = BigInt::inv_mod(pub_e, phi);
        if (priv_d) {
            BigInt msg = BigInt::random_below(mod_n).value();
            BigInt ct = BigInt::pow_mod(msg, pub_e, mod_n).value();
            BigInt pt = BigInt::pow_mod(ct, *priv_d, mod_n).value();
            CHECK(pt == msg);
        }
    }
}

void test_bigint_to_f64() {
    using keel::bignum::BigInt;

    CHECK(BigInt(0).to_f64() == 0.0);
    CHECK(BigInt(-42).to_f64() == -42.0);
    BigInt huge = BigInt(1) << 1025;
    CHECK(keel::is_inf(huge.to_f64()) && huge.to_f64() > 0.0);
    CHECK(keel::is_inf((-huge).to_f64()) && (-huge).to_f64() < 0.0);

    BigInt p200 = BigInt(1) << 200;
    CHECK(p200.to_f64() == keel::ldexp(1.0, 200));
}

void test_rational() {
    using keel::bignum::BigInt;
    using keel::bignum::Rational;

    Rational r0;
    CHECK(r0.is_zero());
    CHECK(r0.is_integer());
    CHECK(r0.to_string() == "0");
    CHECK(!r0.inv().has_value());

    Rational a(BigInt(6), BigInt(-9));
    CHECK(a.num() == BigInt(-2));
    CHECK(a.den() == BigInt(3));
    CHECK(!a.is_integer());
    CHECK(a.to_string() == "-2/3");
    CHECK(a.inv().value() == Rational(-3, 2));

    // floor / ceil / trunc on positive, negative, integer.
    CHECK(Rational(7, 3).floor() == BigInt(2));
    CHECK(Rational(7, 3).ceil() == BigInt(3));
    CHECK(Rational(7, 3).trunc() == BigInt(2));
    CHECK(Rational(-7, 3).floor() == BigInt(-3));
    CHECK(Rational(-7, 3).ceil() == BigInt(-2));
    CHECK(Rational(-7, 3).trunc() == BigInt(-2));
    CHECK(Rational(-6, 3).floor() == BigInt(-2));
    CHECK(Rational(-6, 3).ceil() == BigInt(-2));
    CHECK(Rational(-6, 3).trunc() == BigInt(-2));

    Rational b(1, 3);
    CHECK(a + b == Rational(-1, 3));
    CHECK(a - b == Rational(-1));
    CHECK(a * b == Rational(-2, 9));
    CHECK(a / b == Rational(-2));
    CHECK(Rational::pow_u64(Rational(2, 3), 5) == Rational(32, 243));
    CHECK(Rational::pow_i64(Rational(2, 3), -3).value() == Rational(27, 8));
    CHECK(Rational::pow_i64(Rational(-2, 3), -3).value() == Rational(-27, 8));
    CHECK(Rational::pow_i64(Rational(0), 0).value() == Rational(1));
    CHECK(!Rational::pow_i64(Rational(0), -1).has_value());
    CHECK(Rational::pow_i64(Rational(1), keel::Limits<keel::i64>::kMin).value() == Rational(1));
    CHECK(Rational::pow_i64(Rational(-1), keel::Limits<keel::i64>::kMin).value() == Rational(1));

    keel::i64 out = 0;
    CHECK(Rational(-7, 1).fits_i64(out) && out == -7);
    CHECK(!Rational(-7, 2).fits_i64(out));

    Rational parsed;
    CHECK(Rational::try_parse(parsed, "-12/18") && parsed == Rational(-2, 3));
    CHECK(Rational::try_parse(parsed, "42") && parsed == Rational(42));
    CHECK(!Rational::try_parse(parsed, "1/0"));
    CHECK(!Rational::try_parse(parsed, "1/-2"));
    CHECK(!Rational::try_parse(parsed, "abc"));

    // Decimal / scientific notation parsing.
    CHECK(Rational::try_parse_decimal(parsed, "0.1") && parsed == Rational(1, 10));
    CHECK(Rational::try_parse_decimal(parsed, "-3.25") && parsed == Rational(-13, 4));
    CHECK(Rational::try_parse_decimal(parsed, ".5") && parsed == Rational(1, 2));
    CHECK(Rational::try_parse_decimal(parsed, "5.") && parsed == Rational(5, 1));
    CHECK(Rational::try_parse_decimal(parsed, "1.25e-2") && parsed == Rational(1, 80));
    CHECK(Rational::try_parse_decimal(parsed, "-2.5E+3") && parsed == Rational(-2500, 1));
    CHECK(Rational::try_parse_decimal(parsed, "-0.000") && parsed.is_zero());
    CHECK(!Rational::try_parse_decimal(parsed, ""));
    CHECK(!Rational::try_parse_decimal(parsed, "."));
    CHECK(!Rational::try_parse_decimal(parsed, "1e"));
    CHECK(!Rational::try_parse_decimal(parsed, "1.2.3"));
    CHECK(!Rational::try_parse_decimal(parsed, "1e99999999"));

    // Exact double round-trip including subnormals!
    const keel::f64 samples[] = {
        0.0,
        1.0,
        -3.5,
        0.1,
        1e-100,
        -1e100,
        3.141592653589793,
        bits_to_f64(0x0000000000000001ULL), // DBL_TRUE_MIN = 2^-1074
        bits_to_f64(0x8000000000000001ULL), // -2^-1074
        bits_to_f64(0x000FFFFFFFFFFFFFULL), // largest subnormal
        bits_to_f64(0x0010000000000000ULL), // smallest normal = 2^-1022
        bits_to_f64(0x7FEFFFFFFFFFFFFFULL), // DBL_MAX
    };
    for (keel::f64 x : samples) {
        Rational rx = Rational::from_f64_bits(x);
        CHECK(rx.to_f64() == x);
    }

    // Subnormal tie-to-even single-rounding check:
    // 2^-1075 == half of DBL_TRUE_MIN -> rounds to even (0.0).
    Rational half_min(BigInt(1), BigInt(1) << 1075);
    CHECK(half_min.to_f64() == 0.0);
    // (2^60 + 1) / 2^1135 > 2^-1075 -> rounds UP to DBL_TRUE_MIN (2^-1074).
    Rational just_above_half((BigInt(1) << 60) + BigInt(1), BigInt(1) << 1135);
    CHECK(f64_to_bits(just_above_half.to_f64()) == 0x0000000000000001ULL);

    // Large ratio (num and den both >> 2^1024, ratio == 3/4).
    BigInt big = BigInt(1) << 2000;
    Rational huge_ratio(BigInt(3) * big, BigInt(4) * big);
    CHECK(huge_ratio.to_f64() == 0.75);

    // Huge overflow / underflow fast paths (must not allocate massive shifts).
    Rational massive(BigInt(1) << 5000, BigInt(3));
    CHECK(keel::is_inf(massive.to_f64()) && massive.to_f64() > 0.0);
    Rational tiny(BigInt(1), BigInt(1) << 5000);
    CHECK(tiny.to_f64() == 0.0);
}

// Splits a line on single ASCII spaces into tokens.
keel::Array<keel::StringRef> split_spaces(keel::StringRef line) {
    keel::Array<keel::StringRef> toks;
    keel::usize i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\r' || line[i] == '\t')) ++i;
        if (i >= line.size()) break;
        const keel::usize start = i;
        while (i < line.size() && line[i] != ' ' && line[i] != '\r' && line[i] != '\t') ++i;
        toks.push_back(line.substr(start, i - start));
    }
    return toks;
}

keel::i32 run_batch_file(keel::StringRef path) {
    using keel::bignum::BigInt;
    using keel::bignum::Rational;

    auto res = keel::fs_read_text(path);
    if (res.is_err()) {
        keel::eprintln("failed to read batch file: {}", path);
        return 2;
    }
    const keel::String& text = res.value();
    keel::usize pos = 0;
    while (pos < text.size()) {
        keel::usize end = pos;
        while (end < text.size() && text[end] != '\n') ++end;
        const keel::StringRef line = keel::StringRef(text.data() + pos, end - pos);
        pos = (end < text.size()) ? (end + 1) : end;

        auto toks = split_spaces(line);
        if (toks.empty()) continue;
        const keel::StringRef cmd = toks[0];

        auto parse_bi = [](keel::StringRef s) -> BigInt {
            BigInt v;
            const bool ok = BigInt::try_parse(v, s);
            KEEL_ASSERT(ok);
            return v;
        };
        auto parse_rat = [](keel::StringRef s) -> Rational {
            Rational v;
            const bool ok = Rational::try_parse(v, s);
            KEEL_ASSERT(ok);
            return v;
        };
        auto parse_u64 = [&parse_bi](keel::StringRef s) -> keel::u64 {
            keel::u64 v = 0;
            const bool ok = parse_bi(s).fits_u64(v);
            KEEL_ASSERT(ok);
            return v;
        };
        auto parse_i64 = [&parse_bi](keel::StringRef s) -> keel::i64 {
            keel::i64 v = 0;
            const bool ok = parse_bi(s).fits_i64(v);
            KEEL_ASSERT(ok);
            return v;
        };

        if (cmd == "ADD") {
            keel::println("{}", (parse_bi(toks[1]) + parse_bi(toks[2])).to_string());
        } else if (cmd == "SUB") {
            keel::println("{}", (parse_bi(toks[1]) - parse_bi(toks[2])).to_string());
        } else if (cmd == "MUL") {
            keel::println("{}", (parse_bi(toks[1]) * parse_bi(toks[2])).to_string());
        } else if (cmd == "DIVMOD") {
            BigInt q;
            BigInt r;
            BigInt::div_mod(parse_bi(toks[1]), parse_bi(toks[2]), q, r);
            keel::println("{} {}", q.to_string(), r.to_string());
        } else if (cmd == "MODEUC") {
            keel::println("{}", BigInt::mod_euclid(parse_bi(toks[1]), parse_bi(toks[2])).to_string());
        } else if (cmd == "SHL") {
            keel::println("{}", (parse_bi(toks[1]) << keel::as<keel::usize>(parse_u64(toks[2]))).to_string());
        } else if (cmd == "SHR") {
            keel::println("{}", (parse_bi(toks[1]) >> keel::as<keel::usize>(parse_u64(toks[2]))).to_string());
        } else if (cmd == "AND") {
            keel::println("{}", (parse_bi(toks[1]) & parse_bi(toks[2])).to_string());
        } else if (cmd == "OR") {
            keel::println("{}", (parse_bi(toks[1]) | parse_bi(toks[2])).to_string());
        } else if (cmd == "XOR") {
            keel::println("{}", (parse_bi(toks[1]) ^ parse_bi(toks[2])).to_string());
        } else if (cmd == "NOT") {
            keel::println("{}", (~parse_bi(toks[1])).to_string());
        } else if (cmd == "BITINFO") {
            BigInt a = parse_bi(toks[1]);
            const keel::u64 idx = parse_u64(toks[2]);
            keel::println("{} {} {} {} {}",
                          a.bit_length(),
                          a.byte_length(),
                          a.trailing_zeros(),
                          a.popcount(),
                          a.test_bit(keel::as<keel::usize>(idx)) ? 1 : 0);
        } else if (cmd == "GCD") {
            keel::println("{}", BigInt::gcd(parse_bi(toks[1]), parse_bi(toks[2])).to_string());
        } else if (cmd == "LCM") {
            keel::println("{}", BigInt::lcm(parse_bi(toks[1]), parse_bi(toks[2])).to_string());
        } else if (cmd == "GCDEXT") {
            BigInt x;
            BigInt y;
            BigInt g = BigInt::gcd_ext(parse_bi(toks[1]), parse_bi(toks[2]), x, y);
            keel::println("{} {} {}", g.to_string(), x.to_string(), y.to_string());
        } else if (cmd == "INVMOD") {
            auto r = BigInt::inv_mod(parse_bi(toks[1]), parse_bi(toks[2]));
            if (r) {
                keel::println("{}", r->to_string());
            } else {
                keel::println("NONE");
            }
        } else if (cmd == "POWMOD") {
            auto r = BigInt::pow_mod(parse_bi(toks[1]), parse_bi(toks[2]), parse_bi(toks[3]));
            if (r) {
                keel::println("{}", r->to_string());
            } else {
                keel::println("NONE");
            }
        } else if (cmd == "JACOBI") {
            keel::println("{}", BigInt::jacobi(parse_bi(toks[1]), parse_bi(toks[2])));
        } else if (cmd == "ISQRT") {
            auto r = parse_bi(toks[1]).isqrt();
            if (r) {
                keel::println("{}", r->to_string());
            } else {
                keel::println("NONE");
            }
        } else if (cmd == "RADIX") {
            const keel::i64 in_base = parse_i64(toks[2]);
            const keel::i64 out_base = parse_i64(toks[3]);
            BigInt v;
            if (!BigInt::try_parse_radix(v, toks[1], keel::as<keel::i32>(in_base))) {
                keel::println("ERR");
            } else {
                keel::println("{}", v.to_string_radix(keel::as<keel::i32>(out_base)));
            }
        } else if (cmd == "BYTES_BE") {
            keel::Array<keel::u8> raw;
            hex_to_bytes(toks[1], raw);
            const keel::u64 min_len = parse_u64(toks[2]);
            BigInt v = BigInt::from_bytes_be(keel::ArrayRef<keel::u8>(raw));
            auto be = v.to_bytes_be(keel::as<keel::usize>(min_len));
            auto le = v.to_bytes_le(keel::as<keel::usize>(min_len));
            keel::println("{} {} {}",
                          v.to_string(),
                          be.empty() ? keel::String::from("-") : bytes_to_hex(keel::ArrayRef<keel::u8>(be)),
                          le.empty() ? keel::String::from("-") : bytes_to_hex(keel::ArrayRef<keel::u8>(le)));
        } else if (cmd == "BYTES_LE") {
            keel::Array<keel::u8> raw;
            hex_to_bytes(toks[1], raw);
            const keel::u64 min_len = parse_u64(toks[2]);
            BigInt v = BigInt::from_bytes_le(keel::ArrayRef<keel::u8>(raw));
            auto be = v.to_bytes_be(keel::as<keel::usize>(min_len));
            auto le = v.to_bytes_le(keel::as<keel::usize>(min_len));
            keel::println("{} {} {}",
                          v.to_string(),
                          be.empty() ? keel::String::from("-") : bytes_to_hex(keel::ArrayRef<keel::u8>(be)),
                          le.empty() ? keel::String::from("-") : bytes_to_hex(keel::ArrayRef<keel::u8>(le)));
        } else if (cmd == "BIGF64") {
            keel::println("{}", u64_to_hex16(f64_to_bits(parse_bi(toks[1]).to_f64())));
        } else if (cmd == "RAT_OP") {
            const keel::StringRef op = toks[1];
            const Rational a = parse_rat(toks[2]);
            const Rational b = parse_rat(toks[3]);
            if (op == "+") {
                keel::println("{}", (a + b).to_string());
            } else if (op == "-") {
                keel::println("{}", (a - b).to_string());
            } else if (op == "*") {
                keel::println("{}", (a * b).to_string());
            } else if (op == "/") {
                keel::println("{}", (a / b).to_string());
            } else if (op == "<") {
                keel::println("{}", (a < b) ? 1 : 0);
            }
        } else if (cmd == "RAT_POW") {
            const keel::i64 e = parse_i64(toks[2]);
            auto r = Rational::pow_i64(parse_rat(toks[1]), e);
            if (r) {
                keel::println("{}", r->to_string());
            } else {
                keel::println("NONE");
            }
        } else if (cmd == "RAT_ROUND") {
            const Rational a = parse_rat(toks[1]);
            keel::println("{} {} {}", a.floor().to_string(), a.ceil().to_string(), a.trunc().to_string());
        } else if (cmd == "RAT_F64") {
            keel::println("{}", u64_to_hex16(f64_to_bits(parse_rat(toks[1]).to_f64())));
        } else if (cmd == "RAT_FROM_F64") {
            keel::u64 bits = 0;
            parse_u64_hex(toks[1], bits);
            keel::println("{}", Rational::from_f64_bits(bits_to_f64(bits)).to_string());
        } else if (cmd == "RAT_DEC") {
            Rational r;
            if (!Rational::try_parse_decimal(r, toks[1])) {
                keel::println("ERR");
            } else {
                keel::println("{}", r.to_string());
            }
        } else if (cmd == "PRIME") {
            keel::println("{}", BigInt::is_prime(parse_bi(toks[1])) ? 1 : 0);
        } else if (cmd == "FACTORIZE") {
            auto fac = BigInt::factorize(parse_bi(toks[1]));
            keel::String line_out;
            for (keel::usize i = 0; i < fac.size(); ++i) {
                if (i > 0) line_out.push_back(' ');
                line_out.append(fac[i].prime.to_string());
                line_out.push_back('^');
                char ebuf[16];
                const keel::usize en = keel::format_i64(ebuf, sizeof(ebuf), keel::as<keel::i64>(fac[i].exponent));
                line_out.append(keel::StringRef(ebuf, en));
            }
            keel::println("{}", line_out.empty() ? keel::StringRef("1") : line_out.as_ref());
        } else if (cmd == "GENPRIME") {
            const keel::u64 bits = parse_u64(toks[1]);
            auto p = BigInt::random_prime(keel::as<keel::usize>(bits));
            if (p) {
                keel::println("{}", p->to_string());
            } else {
                keel::println("NONE");
            }
        } else {
            keel::eprintln("unknown batch command: {}", cmd);
            return 2;
        }
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 3 && keel::StringRef(argv[1]) == "--batch") {
        return run_batch_file(keel::StringRef(argv[2]));
    }

    test_bigint_basic_and_sbo();
    test_bigint_arith_div_shift_karatsuba();
    test_bigint_bitwise();
    test_bigint_radix_hex_and_bytes();
    test_bigint_modular_and_number_theory();
    test_bigint_primality_factorize_and_csprng();
    test_bigint_to_f64();
    test_rational();

    if (g_failed != 0) {
        keel::eprintln("test_bignum: {} failure(s)", g_failed);
        return 1;
    }
    keel::println("test_bignum: all tests passed");
    return 0;
}
