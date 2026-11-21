#include <keel/bignum.h>
#include <keel/math.h>

#include <cstdio>

namespace {

keel::i32 g_failed = 0;

void check(bool cond, const char* expr, keel::i32 line) {
    if (!cond) {
        std::fprintf(stderr, "FAIL [line %d]: %s\n", static_cast<int>(line), expr);
        ++g_failed;
    }
}

#define CHECK(cond) check((cond), #cond, __LINE__)

void test_bigint_basic_and_sbo() {
    using keel::BigInt;

    BigInt zero;
    CHECK(zero.is_zero());
    CHECK(zero.is_inline());
    CHECK(zero.signum() == 0);
    CHECK(zero.bit_length() == 0);
    CHECK(zero.to_string() == "0");

    BigInt one(1);
    BigInt neg_min(keel::Limits<keel::i64>::kMin);
    CHECK(one.is_inline());
    CHECK(neg_min.is_inline());
    CHECK(neg_min.is_negative());
    CHECK(neg_min.to_string() == "-9223372036854775808");

    keel::i64 i_out = 0;
    CHECK(neg_min.fits_i64(i_out) && i_out == keel::Limits<keel::i64>::kMin);
    keel::u64 u_out = 0;
    CHECK(!neg_min.fits_u64(u_out));

    // 2 limbs inline, 3 limbs heap.
    BigInt two_limbs = (BigInt(1) << 127) - BigInt(1);
    CHECK(two_limbs.is_inline());
    CHECK(two_limbs.bit_length() == 127);

    BigInt three_limbs = BigInt(1) << 128;
    CHECK(!three_limbs.is_inline());
    CHECK(three_limbs.bit_length() == 129);

    BigInt shrunk = three_limbs >> 64;
    BigInt copy_shrunk = shrunk;
    CHECK(copy_shrunk.is_inline());
    CHECK(copy_shrunk == (BigInt(1) << 64));
}

void test_bigint_arith_div_shift() {
    using keel::BigInt;

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

    BigInt::div_mod(-a, b, q, r);
    CHECK(q * b + r == -a);
    CHECK(r.is_negative() || r.is_zero());

    // Knuth D borrow edge case: high-bit limbs.
    BigInt u = (BigInt(1) << 192) - BigInt(1);
    BigInt v = (BigInt(1) << 96) + BigInt(1);
    BigInt::div_mod(u, v, q, r);
    CHECK(q * v + r == u);
    CHECK(q == (BigInt(1) << 96) - BigInt(1));
    CHECK(r.is_zero());

    // Negative right shift uses floor semantics.
    CHECK((BigInt(-1) >> 1) == BigInt(-1));
    CHECK((BigInt(-5) >> 1) == BigInt(-3));
    CHECK((BigInt(-4) >> 1) == BigInt(-2));

    // GCD / LCM / pow / sqrt_exact.
    CHECK(BigInt::gcd(BigInt(48), BigInt(-18)) == BigInt(6));
    CHECK(BigInt::lcm(BigInt(12), BigInt(18)) == BigInt(36));
    CHECK(BigInt::pow_u64(BigInt(0), 0) == BigInt(1));
    CHECK(BigInt::pow_u64(BigInt(3), 40).to_string() == "12157665459056928801");

    BigInt sq_root;
    BigInt sq = a * a;
    CHECK(sq.sqrt_exact(sq_root) && sq_root == a);
    CHECK(!(sq + BigInt(1)).sqrt_exact(sq_root));
    CHECK(!BigInt(-4).sqrt_exact(sq_root));

    // Invalid parse inputs.
    BigInt bad;
    CHECK(!BigInt::try_parse(bad, ""));
    CHECK(!BigInt::try_parse(bad, "-"));
    CHECK(!BigInt::try_parse(bad, " 123"));
    CHECK(!BigInt::try_parse(bad, "12a"));
}

void test_bigint_primality_and_factorize() {
    using keel::BigInt;

    CHECK(!BigInt::is_prime(BigInt(-7)));
    CHECK(!BigInt::is_prime(BigInt(0)));
    CHECK(!BigInt::is_prime(BigInt(1)));
    CHECK(BigInt::is_prime(BigInt(2)));
    CHECK(BigInt::is_prime(BigInt(3)));
    CHECK(!BigInt::is_prime(BigInt(561)));       // Carmichael number
    CHECK(!BigInt::is_prime(BigInt(341550071728321ULL))); // strong pseudoprime to 2..17
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
        prod *= BigInt::pow_u64(fac[i].prime, static_cast<keel::u64>(fac[i].exponent));
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
}

void test_bigint_to_f64() {
    using keel::BigInt;

    CHECK(BigInt(0).to_f64() == 0.0);
    CHECK(BigInt(-42).to_f64() == -42.0);
    // Overflow past 2^1024 -> +/-inf.
    BigInt huge = BigInt(1) << 1025;
    CHECK(keel::is_inf(huge.to_f64()) && huge.to_f64() > 0.0);
    CHECK(keel::is_inf((-huge).to_f64()) && (-huge).to_f64() < 0.0);

    // Exact power of 2 near 2^200.
    BigInt p200 = BigInt(1) << 200;
    CHECK(p200.to_f64() == keel::ldexp(1.0, 200));
}

void test_rational() {
    using keel::BigInt;
    using keel::Rational;

    Rational r0;
    CHECK(r0.is_zero());
    CHECK(r0.to_string() == "0");

    Rational a(BigInt(6), BigInt(-9));
    CHECK(a.num() == BigInt(-2));
    CHECK(a.den() == BigInt(3));
    CHECK(a.to_string() == "-2/3");

    Rational b(1, 3);
    CHECK(a + b == Rational(-1, 3));
    CHECK(a - b == Rational(-1));
    CHECK(a * b == Rational(-2, 9));
    CHECK(a / b == Rational(-2));
    CHECK(Rational::pow_u64(Rational(2, 3), 5) == Rational(32, 243));

    keel::i64 out = 0;
    CHECK(Rational(-7, 1).fits_i64(out) && out == -7);
    CHECK(!Rational(-7, 2).fits_i64(out));

    Rational parsed;
    CHECK(Rational::try_parse(parsed, "-12/18") && parsed == Rational(-2, 3));
    CHECK(Rational::try_parse(parsed, "42") && parsed == Rational(42));
    CHECK(!Rational::try_parse(parsed, "1/0"));
    CHECK(!Rational::try_parse(parsed, "1/-2"));
    CHECK(!Rational::try_parse(parsed, "abc"));

    // Exact double round-trip and scale-invariant to_f64.
    const keel::f64 samples[] = {0.0, 1.0, -3.5, 0.1, 1e-100, -1e100, 3.141592653589793};
    for (keel::f64 x : samples) {
        Rational rx = Rational::from_f64_bits(x);
        CHECK(rx.to_f64() == x);
    }

    // Large ratio (num and den both >> 2^1024, ratio == 3/4).
    BigInt big = BigInt(1) << 2000;
    Rational huge_ratio(BigInt(3) * big, BigInt(4) * big);
    CHECK(huge_ratio.to_f64() == 0.75);

    // Non-reducible large ratio near 1/3.
    BigInt big_num = (BigInt(1) << 1500) + BigInt(1);
    BigInt big_den = BigInt(3) * (BigInt(1) << 1500) + BigInt(2);
    Rational near_third(big_num, big_den);
    CHECK(keel::abs(near_third.to_f64() - (1.0 / 3.0)) < 1e-15);
}

} // namespace

int main() {
    test_bigint_basic_and_sbo();
    test_bigint_arith_div_shift();
    test_bigint_primality_and_factorize();
    test_bigint_to_f64();
    test_rational();

    if (g_failed != 0) {
        std::fprintf(stderr, "test_bignum: %d failure(s)\n", static_cast<int>(g_failed));
        return 1;
    }
    std::printf("test_bignum: all tests passed\n");
    return 0;
}
