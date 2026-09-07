#include "science.hpp"
#include "performance_diagnostics.hpp"

#include <amsolcpp/coordinates.hpp>
#include <amsolcpp/exceptions.hpp>

#include <cstdint>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace amsolcpp::detail {
namespace {

constexpr double bohr_radius_angstrom = 0.529167;
constexpr double ev_per_atomic_unit = 27.21;
constexpr double rt3 = 0.57735026918962576451;
constexpr std::array<std::pair<int, int>, 10> product_pairs{{
    {0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1},
    {2, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}}};
constexpr std::array<int, 9> l1scat{{2, 4, 7, 3, 5, 6, 8, 9, 10}};
constexpr std::array<int, 99> l2scat{{
    11,31,61,21,41,51,71,81,91, 2,4,7,3,5,6,8,9,10,
    12,32,34,62,64,67,14,17,37, 22,42,52,72,82,92,13,15,16,18,19,20,
    24,44,54,74,84,94,33,35,36,38,39,40,
    27,47,57,77,87,97,63,65,66,68,69,70,
    23,25,26,28,29,30, 43,45,46,48,49,50, 53,55,56,58,59,60,
    73,75,76,78,79,80, 83,85,86,88,89,90, 93,95,96,98,99,100}};

double distance(const Atom& a, const Atom& b) {
    return std::hypot(a.x_angstrom - b.x_angstrom,
                      a.y_angstrom - b.y_angstrom,
                      a.z_angstrom - b.z_angstrom);
}

int principal_shell(const int z) {
    if (z == 1) return 1;
    if (z <= 9) return 2;
    if (z <= 17) return 3;
    if (z == 35) return 4;
    return 5;
}

double coefficient(const int index) {
    return index % 2 == 0 ? 0.0 : 2.0 / static_cast<double>(index);
}

struct SetValues {
    double sa{};
    double sb{};
    int isp{};
    int ips{};
    std::array<double, 8> a{};
    std::array<double, 8> b{};
};

SetValues set_integrals(double s1, double s2, const int na, const int nb,
                        const double rab, const int ii) {
    SetValues result;
    if (na <= nb) {
        result.isp = 1; result.ips = 2; result.sa = s1; result.sb = s2;
    } else {
        result.isp = 2; result.ips = 1; result.sa = s2; result.sb = s1;
    }
    const int k = ii > 3 ? ii : ii + 1;
    double x = 0.5 * rab * (result.sa + result.sb);
    const double y = 1.0 / x;
    const double ex = std::exp(-x);
    result.a[1] = ex * y;
    for (int index = 1; index <= k; ++index) {
        result.a[static_cast<std::size_t>(index + 1)] =
            (result.a[static_cast<std::size_t>(index)] * static_cast<double>(index) + ex) * y;
    }
    x = 0.5 * rab * (result.sb - result.sa);
    if (std::abs(x) > 3.0) {
        double exp_x = std::exp(x);
        const double exp_minus_x = 1.0 / exp_x;
        const double inverse_x = 1.0 / x;
        result.b[1] = (exp_x - exp_minus_x) * inverse_x;
        for (int index = 1; index <= k; ++index) {
            exp_x = -exp_x;
            result.b[static_cast<std::size_t>(index + 1)] =
                (static_cast<double>(index) * result.b[static_cast<std::size_t>(index)] +
                 exp_x - exp_minus_x) * inverse_x;
        }
    } else {
        for (int index = 1; index <= k + 1; ++index) result.b[static_cast<std::size_t>(index)] = coefficient(index);
        if (std::abs(x) > 1.0e-6) {
            const int last = static_cast<int>(5.0 + 4.0 * std::abs(x));
            std::vector<double> work(static_cast<std::size_t>(last + 1), 0.0);
            double inverse_factorial = 1.0;
            for (int m = 1; m <= last; ++m) {
                inverse_factorial /= static_cast<double>(m);
                work[static_cast<std::size_t>(m)] = inverse_factorial * std::pow(-x, m);
            }
            for (int index = 1; index <= k + 1; ++index) {
                for (int m = index % 2 + 1; m <= last; m += 2) {
                    result.b[static_cast<std::size_t>(index)] +=
                        work[static_cast<std::size_t>(m)] * coefficient(m + index);
                }
            }
        }
    }
    return result;
}

std::array<double, 14> bfn(const double x) {
    std::array<double, 14> values{};
    if (std::abs(x) > 3.0) {
        double exp_x = std::exp(x);
        const double exp_minus_x = 1.0 / exp_x;
        const double inverse_x = 1.0 / x;
        values[1] = (exp_x - exp_minus_x) * inverse_x;
        for (int index = 1; index < 13; ++index) {
            exp_x = -exp_x;
            values[static_cast<std::size_t>(index + 1)] =
                (static_cast<double>(index) * values[static_cast<std::size_t>(index)] +
                 exp_x - exp_minus_x) * inverse_x;
        }
        return values;
    }
    for (int index = 1; index < 14; ++index) values[static_cast<std::size_t>(index)] = coefficient(index);
    if (std::abs(x) <= 1.0e-6) return values;
    const int last = static_cast<int>(5.0 + 4.0 * std::abs(x));
    double inverse_factorial = 1.0;
    for (int m = 1; m <= last; ++m) {
        inverse_factorial /= static_cast<double>(m);
        const double work = std::pow(-x, m) * inverse_factorial;
        for (int index = 1; index < 14; ++index) {
            if (m >= index % 2 + 1 && (m - (index % 2 + 1)) % 2 == 0) {
                values[static_cast<std::size_t>(index)] += work * coefficient(m + index);
            }
        }
    }
    return values;
}

double combination(const int n, const int k) {
    if (k < 0 || k > n) return 0.0;
    double value = 1.0;
    for (int i = 1; i <= std::min(k, n - k); ++i) {
        value *= static_cast<double>(n - i + 1) / static_cast<double>(i);
    }
    return value;
}

double ss_overlap(const int na, const int nb, const int la, const int lb,
                  const int m, const double ua, const double ub, const double dist) {
    std::array<double, 28> factorial{};
    factorial[0] = 1.0;
    for (std::size_t i = 1; i < factorial.size(); ++i) factorial[i] = factorial[i - 1] * static_cast<double>(i);
    const double rab = dist / bohr_radius_angstrom;
    const double p = (ua + ub) * rab * 0.5;
    const double ba = (ua - ub) * rab * 0.5;
    std::array<double, 21> af{};
    af[1] = std::exp(-p) / p;
    for (int i = 1; i < 20; ++i) af[static_cast<std::size_t>(i + 1)] = static_cast<double>(i) * af[static_cast<std::size_t>(i)] / p + af[1];
    const auto bf = bfn(ba);
    auto angular = [](const int l, const int mm, const int i) {
        if (l == 1 && mm == 1 && i == 1) return 1.0;
        if (l == 2 && mm == 1 && i == 1) return 1.0;
        if (l == 2 && mm == 2 && i == 1) return 1.0;
        if (l == 3 && mm == 1 && i == 1) return 1.5;
        if (l == 3 && mm == 2 && i == 1) return 1.73205;
        if (l == 3 && mm == 3 && i == 1) return 1.224745;
        if (l == 3 && mm == 1 && i == 3) return -0.5;
        return 0.0;
    };
    double total = 0.0;
    for (int iv = 1; iv <= la - m + 1; iv += 2) {
        const int av = na + iv - la;
        const int ic = la + 2 - iv - m;
        for (int jv = 1; jv <= lb - m + 1; jv += 2) {
            const int bv = nb + jv - lb;
            const int id = lb - jv - m + 2;
            double subtotal = 0.0;
            const int ia = av + 1, ib = bv + 1, ab = av + bv - 1;
            for (int k1 = 1; k1 <= ia; ++k1) for (int k2 = 1; k2 <= ib; ++k2)
            for (int k3 = 1; k3 <= ic; ++k3) for (int k4 = 1; k4 <= id; ++k4)
            for (int k5 = 1; k5 <= m; ++k5) for (int k6 = 1; k6 <= m; ++k6) {
                const double part = combination(ia - 1, k1 - 1) * combination(ib - 1, k2 - 1) *
                    combination(ic - 1, k3 - 1) * combination(id - 1, k4 - 1) *
                    combination(m - 1, k5 - 1) * combination(m - 1, k6 - 1);
                const int q = ab - k1 - k2 + k3 + k4 + 2 * k5;
                const int pp = k1 + k2 + k3 + k4 + 2 * k6 - 5;
                const int jx = m + k2 + k4 + k5 + k6 - 5;
                const int ix = jx / 2;
                subtotal += part * (static_cast<double>(ix * 2 - jx) + 0.5) *
                    af[static_cast<std::size_t>(q)] * bf[static_cast<std::size_t>(pp)];
            }
            total += subtotal * angular(la, m, iv) * angular(lb, m, jv) * 2.0;
        }
    }
    const double scale = std::pow(rab, na + nb + 1) * std::pow(ua, na) * std::pow(ub, nb);
    const double norm = std::sqrt(ua * ub / (factorial[static_cast<std::size_t>(2 * na)] *
        factorial[static_cast<std::size_t>(2 * nb)]) * static_cast<double>((2 * la - 1) * (2 * lb - 1)));
    return total * scale * norm / std::pow(2.0, m);
}

using Primitives = std::array<std::array<std::array<double, 4>, 3>, 3>;

Primitives primitives(const ElementParameters& a, const ElementParameters& b, const double dist) {
    Primitives s{};
    if (dist < 0.1 || dist >= 10.0) return s;
    const int shell_a = principal_shell(a.atomic_number), shell_b = principal_shell(b.atomic_number);
    if (shell_a > 2 || shell_b > 2) {
        const std::array<double, 3> ea{{0.0, a.zs, a.zp}};
        const std::array<double, 3> eb{{0.0, b.zs, b.zp}};
        s[1][1][1] = ss_overlap(shell_a, shell_b, 1, 1, 1, ea[1], eb[1], dist);
        s[1][2][1] = ss_overlap(shell_a, shell_b, 1, 2, 1, ea[1], eb[2], dist);
        s[2][1][1] = ss_overlap(shell_a, shell_b, 2, 1, 1, ea[2], eb[1], dist);
        s[2][2][1] = ss_overlap(shell_a, shell_b, 2, 2, 1, ea[2], eb[2], dist);
        s[2][2][2] = ss_overlap(shell_a, shell_b, 2, 2, 2, ea[2], eb[2], dist);
        return s;
    }
    const double rab = dist / bohr_radius_angstrom;
    if (a.atomic_number == 1 && b.atomic_number == 1) {
        const auto q = set_integrals(a.zs, b.zs, a.atomic_number, b.atomic_number, rab, 1);
        s[1][1][1] = 0.25 * std::sqrt(std::pow(q.sa * q.sb * rab * rab, 3)) *
            (q.a[3] * q.b[1] - q.b[3] * q.a[1]);
        return s;
    }
    if (a.atomic_number == 1 || b.atomic_number == 1) {
        auto q = set_integrals(a.zs, b.zs, a.atomic_number, b.atomic_number, rab, 2);
        double w = std::sqrt(std::pow(q.sa, 3) * std::pow(q.sb, 5)) * std::pow(rab, 4) * 0.125;
        s[1][1][1] = w * rt3 * (q.a[4] * q.b[1] - q.b[4] * q.a[1] + q.a[3] * q.b[2] - q.b[3] * q.a[2]);
        q = b.atomic_number > 1
            ? set_integrals(a.zs, b.zp, a.atomic_number, b.atomic_number, rab, 2)
            : set_integrals(a.zp, b.zs, a.atomic_number, b.atomic_number, rab, 2);
        w = std::sqrt(std::pow(q.sa, 3) * std::pow(q.sb, 5)) * std::pow(rab, 4) * 0.125;
        s[static_cast<std::size_t>(q.isp)][static_cast<std::size_t>(q.ips)][1] =
            w * (q.a[3] * q.b[1] - q.b[3] * q.a[1] + q.a[4] * q.b[2] - q.b[4] * q.a[2]);
        return s;
    }
    auto q = set_integrals(a.zs, b.zs, a.atomic_number, b.atomic_number, rab, 4);
    double w = std::sqrt(std::pow(q.sa * q.sb, 5)) * std::pow(rab, 5) * 0.0625;
    s[1][1][1] = w * (q.a[5] * q.b[1] + q.b[5] * q.a[1] - 2.0 * q.a[3] * q.b[3]) / 3.0;
    q = a.atomic_number <= b.atomic_number
        ? set_integrals(a.zs, b.zp, a.atomic_number, b.atomic_number, rab, 4)
        : set_integrals(a.zp, b.zs, a.atomic_number, b.atomic_number, rab, 4);
    w = std::sqrt(std::pow(q.sa * q.sb, 5)) * std::pow(rab, 5) * 0.0625;
    double d = q.a[4] * (q.b[1] - q.b[3]) - q.a[2] * (q.b[3] - q.b[5]);
    double e = q.b[4] * (q.a[1] - q.a[3]) - q.b[2] * (q.a[3] - q.a[5]);
    s[static_cast<std::size_t>(q.isp)][static_cast<std::size_t>(q.ips)][1] = w * rt3 * (d + e);
    q = a.atomic_number <= b.atomic_number
        ? set_integrals(a.zp, b.zs, a.atomic_number, b.atomic_number, rab, 4)
        : set_integrals(a.zs, b.zp, a.atomic_number, b.atomic_number, rab, 4);
    w = std::sqrt(std::pow(q.sa * q.sb, 5)) * std::pow(rab, 5) * 0.0625;
    d = q.a[4] * (q.b[1] - q.b[3]) - q.a[2] * (q.b[3] - q.b[5]);
    e = q.b[4] * (q.a[1] - q.a[3]) - q.b[2] * (q.a[3] - q.a[5]);
    s[static_cast<std::size_t>(q.ips)][static_cast<std::size_t>(q.isp)][1] = -w * rt3 * (e - d);
    q = set_integrals(a.zp, b.zp, a.atomic_number, b.atomic_number, rab, 4);
    w = std::sqrt(std::pow(q.sa * q.sb, 5)) * std::pow(rab, 5) * 0.0625;
    s[2][2][1] = -w * (q.b[3] * (q.a[5] + q.a[1]) - q.a[3] * (q.b[5] + q.b[1]));
    s[2][2][2] = 0.5 * w * (q.a[5] * (q.b[1] - q.b[3]) - q.b[5] * (q.a[1] - q.a[3]) - q.a[3] * q.b[1] + q.b[3] * q.a[1]);
    return s;
}

std::vector<double> overlap_block(const Atom& atom_a, const Atom& atom_b,
                                  const ElementParameters& a, const ElementParameters& b) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::OneElectronTerms);
    const double dx = atom_b.x_angstrom - atom_a.x_angstrom;
    const double dy = atom_b.y_angstrom - atom_a.y_angstrom;
    const double dz = atom_b.z_angstrom - atom_a.z_angstrom;
    const double xy = std::hypot(dx, dy);
    const double r = std::hypot(xy, dz);
    const auto s = primitives(a, b, r);
    double ca{}, cb{}, sa{}, sb{};
    if (xy < 1.0e-10) {
        ca = dz < 0.0 ? -1.0 : (dz > 0.0 ? 1.0 : 0.0); cb = ca;
    } else {
        ca = dx / xy; cb = dz / r; sa = dy / xy; sb = xy / r;
    }
    std::array<double, 76> c{};
    c[37] = 1.0; c[56] = ca * cb; c[41] = ca * sb; c[26] = -sa;
    c[53] = -sb; c[38] = cb; c[23] = 0.0; c[50] = sa * cb; c[35] = sa * sb; c[20] = ca;
    const auto ci = [&c](const int plane, const int first, const int second) {
        return c[static_cast<std::size_t>(first + (second - 1) * 3 + (plane - 1) * 15)];
    };
    constexpr std::array<int, 15> ival{{1,0,9, 1,3,8, 1,4,7, 1,2,6, 0,0,5}};
    std::vector<double> block(static_cast<std::size_t>(a.basis_orbitals * b.basis_orbitals), 0.0);
    const int pqa = a.atomic_number == 1 ? 1 : std::min(principal_shell(a.atomic_number), 3);
    const int pqb = b.atomic_number == 1 ? 1 : std::min(principal_shell(b.atomic_number), 3);
    for (int ii = 1; ii <= pqa; ++ii) for (int jj = 1; jj <= pqb; ++jj) {
        const int aa = jj == 2 ? -1 : 1;
        const int bb = jj == 3 ? -1 : 1;
        for (int k = 4 - ii; k <= 2 + ii; ++k) for (int l = 4 - jj; l <= 2 + jj; ++l) {
            const int row = ival[static_cast<std::size_t>((k - 1) * 3 + ii - 1)];
            const int col = ival[static_cast<std::size_t>((l - 1) * 3 + jj - 1)];
            if (row < 1 || col < 1 || row > a.basis_orbitals || col > b.basis_orbitals) continue;
            const double value = s[static_cast<std::size_t>(ii)][static_cast<std::size_t>(jj)][1] * ci(3, ii, k) * ci(3, jj, l) * aa +
                (ci(4, ii, k) * ci(4, jj, l) + ci(2, ii, k) * ci(2, jj, l)) * bb * s[static_cast<std::size_t>(ii)][static_cast<std::size_t>(jj)][2] +
                (ci(5, ii, k) * ci(5, jj, l) + ci(1, ii, k) * ci(1, jj, l)) * s[static_cast<std::size_t>(ii)][static_cast<std::size_t>(jj)][3];
            block[static_cast<std::size_t>((row - 1) * b.basis_orbitals + col - 1)] = value;
        }
    }
    return block;
}

struct ReppValues {
    std::array<double, 23> ri{};
    std::array<std::array<double, 3>, 5> core{};
};

double invsqrt(const double x) { return 1.0 / std::sqrt(x); }

ReppValues repp(const double r, const ElementParameters& a, const ElementParameters& b) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::TwoElectronIntegrals);
    ReppValues v;
    const double da=a.dd, db=b.dd, qa=a.qq, qb=b.qq;
    const double ee = invsqrt(r*r + 0.25*std::pow(1.0/a.am + 1.0/b.am, 2));
    v.ri[1]=ee*ev_per_atomic_unit;
    v.core[1][1]=-static_cast<double>(b.core_charge)*v.ri[1];
    v.core[1][2]=-static_cast<double>(a.core_charge)*v.ri[1];
    if (a.basis_orbitals < 3 && b.basis_orbitals < 3) return v;
    if (a.basis_orbitals >= 3) {
        const double ade=.25*std::pow(1.0/a.ad+1.0/b.am,2), aqe=.25*std::pow(1.0/a.aq+1.0/b.am,2);
        const double dze=(-invsqrt((r+da)*(r+da)+ade)+invsqrt((r-da)*(r-da)+ade))*.5;
        const double qzze=(invsqrt((r-2*qa)*(r-2*qa)+aqe)-2*invsqrt(r*r+aqe)+invsqrt((r+2*qa)*(r+2*qa)+aqe))*.25;
        const double qxxe=(2*invsqrt(r*r+4*qa*qa+aqe)-2*invsqrt(r*r+aqe))*.25;
        v.ri[2]=-dze; v.ri[3]=ee+qzze; v.ri[4]=ee+qxxe;
        if (b.basis_orbitals < 3) goto finish;
    }
    {
        const double aed=.25*std::pow(1.0/a.am+1.0/b.ad,2), aeq=.25*std::pow(1.0/a.am+1.0/b.aq,2);
        const double edz=(-invsqrt((r-db)*(r-db)+aed)+invsqrt((r+db)*(r+db)+aed))*.5;
        const double eqzz=(invsqrt((r-2*qb)*(r-2*qb)+aeq)-2*invsqrt(r*r+aeq)+invsqrt((r+2*qb)*(r+2*qb)+aeq))*.25;
        const double eqxx=(2*invsqrt(r*r+4*qb*qb+aeq)-2*invsqrt(r*r+aeq))*.25;
        v.ri[5]=-edz; v.ri[11]=ee+eqzz; v.ri[12]=ee+eqxx;
        if (a.basis_orbitals < 3) goto finish;
        const double add=.25*std::pow(1.0/a.ad+1.0/b.ad,2), adq=.25*std::pow(1.0/a.ad+1.0/b.aq,2);
        const double aqd=.25*std::pow(1.0/a.aq+1.0/b.ad,2), aqq=.25*std::pow(1.0/a.aq+1.0/b.aq,2);
        double dxdx=(2*invsqrt(r*r+(da-db)*(da-db)+add)-2*invsqrt(r*r+(da+db)*(da+db)+add))*.25;
        double dzdz=(invsqrt((r+da-db)*(r+da-db)+add)+invsqrt((r-da+db)*(r-da+db)+add)-invsqrt((r-da-db)*(r-da-db)+add)-invsqrt((r+da+db)*(r+da+db)+add))*.25;
        double dzqxx=(-2*invsqrt((r+da)*(r+da)+4*qb*qb+adq)+2*invsqrt((r-da)*(r-da)+4*qb*qb+adq)+2*invsqrt((r+da)*(r+da)+adq)-2*invsqrt((r-da)*(r-da)+adq))*.125;
        double qxxdz=(-2*invsqrt((r-db)*(r-db)+4*qa*qa+aqd)+2*invsqrt((r+db)*(r+db)+4*qa*qa+aqd)+2*invsqrt((r-db)*(r-db)+aqd)-2*invsqrt((r+db)*(r+db)+aqd))*.125;
        double dzqzz=(-invsqrt((r+da-2*qb)*(r+da-2*qb)+adq)+invsqrt((r-da-2*qb)*(r-da-2*qb)+adq)-invsqrt((r+da+2*qb)*(r+da+2*qb)+adq)+invsqrt((r-da+2*qb)*(r-da+2*qb)+adq)-2*invsqrt((r-da)*(r-da)+adq)+2*invsqrt((r+da)*(r+da)+adq))*.125;
        double qzzdz=(-invsqrt((r+2*qa-db)*(r+2*qa-db)+aqd)+invsqrt((r+2*qa+db)*(r+2*qa+db)+aqd)-invsqrt((r-2*qa-db)*(r-2*qa-db)+aqd)+invsqrt((r-2*qa+db)*(r-2*qa+db)+aqd)+2*invsqrt((r-db)*(r-db)+aqd)-2*invsqrt((r+db)*(r+db)+aqd))*.125;
        double qxxqxx=(2*invsqrt(r*r+4*(qa-qb)*(qa-qb)+aqq)+2*invsqrt(r*r+4*(qa+qb)*(qa+qb)+aqq)-4*invsqrt(r*r+4*qa*qa+aqq)-4*invsqrt(r*r+4*qb*qb+aqq)+4*invsqrt(r*r+aqq))*.0625;
        double qxxqyy=(4*invsqrt(r*r+4*qa*qa+4*qb*qb+aqq)-4*invsqrt(r*r+4*qa*qa+aqq)-4*invsqrt(r*r+4*qb*qb+aqq)+4*invsqrt(r*r+aqq))*.0625;
        double qxxqzz=(2*invsqrt((r-2*qb)*(r-2*qb)+4*qa*qa+aqq)+2*invsqrt((r+2*qb)*(r+2*qb)+4*qa*qa+aqq)-2*invsqrt((r-2*qb)*(r-2*qb)+aqq)-2*invsqrt((r+2*qb)*(r+2*qb)+aqq)-4*invsqrt(r*r+4*qa*qa+aqq)+4*invsqrt(r*r+aqq))*.0625;
        double qzzqxx=(2*invsqrt((r+2*qa)*(r+2*qa)+4*qb*qb+aqq)+2*invsqrt((r-2*qa)*(r-2*qa)+4*qb*qb+aqq)-2*invsqrt((r+2*qa)*(r+2*qa)+aqq)-2*invsqrt((r-2*qa)*(r-2*qa)+aqq)-4*invsqrt(r*r+4*qb*qb+aqq)+4*invsqrt(r*r+aqq))*.0625;
        double qzzqzz=(invsqrt((r+2*qa-2*qb)*(r+2*qa-2*qb)+aqq)+invsqrt((r+2*qa+2*qb)*(r+2*qa+2*qb)+aqq)+invsqrt((r-2*qa-2*qb)*(r-2*qa-2*qb)+aqq)+invsqrt((r-2*qa+2*qb)*(r-2*qa+2*qb)+aqq)-2*invsqrt((r-2*qa)*(r-2*qa)+aqq)-2*invsqrt((r+2*qa)*(r+2*qa)+aqq)-2*invsqrt((r-2*qb)*(r-2*qb)+aqq)-2*invsqrt((r+2*qb)*(r+2*qb)+aqq)+4*invsqrt(r*r+aqq))*.0625;
        double dxqxz=(-2*invsqrt((r-qb)*(r-qb)+(da-qb)*(da-qb)+adq)+2*invsqrt((r+qb)*(r+qb)+(da-qb)*(da-qb)+adq)+2*invsqrt((r-qb)*(r-qb)+(da+qb)*(da+qb)+adq)-2*invsqrt((r+qb)*(r+qb)+(da+qb)*(da+qb)+adq))*.125;
        double qxzdx=(-2*invsqrt((r+qa)*(r+qa)+(qa-db)*(qa-db)+aqd)+2*invsqrt((r-qa)*(r-qa)+(qa-db)*(qa-db)+aqd)+2*invsqrt((r+qa)*(r+qa)+(qa+db)*(qa+db)+aqd)-2*invsqrt((r-qa)*(r-qa)+(qa+db)*(qa+db)+aqd))*.125;
        [[maybe_unused]] const double qxyqxy=(4*invsqrt(r*r+2*(qa-qb)*(qa-qb)+aqq)+4*invsqrt(r*r+2*(qa+qb)*(qa+qb)+aqq)-8*invsqrt(r*r+2*(qa*qa+qb*qb)+aqq))*.0625;
        double qxzqxz=(2*invsqrt((r+qa-qb)*(r+qa-qb)+(qa-qb)*(qa-qb)+aqq)-2*invsqrt((r+qa+qb)*(r+qa+qb)+(qa-qb)*(qa-qb)+aqq)-2*invsqrt((r-qa-qb)*(r-qa-qb)+(qa-qb)*(qa-qb)+aqq)+2*invsqrt((r-qa+qb)*(r-qa+qb)+(qa-qb)*(qa-qb)+aqq)-2*invsqrt((r+qa-qb)*(r+qa-qb)+(qa+qb)*(qa+qb)+aqq)+2*invsqrt((r+qa+qb)*(r+qa+qb)+(qa+qb)*(qa+qb)+aqq)+2*invsqrt((r-qa-qb)*(r-qa-qb)+(qa+qb)*(qa+qb)+aqq)-2*invsqrt((r-qa+qb)*(r-qa+qb)+(qa+qb)*(qa+qb)+aqq))*.0625;
        v.ri[6]=dzdz; v.ri[7]=dxdx; v.ri[8]=-edz-qzzdz; v.ri[9]=-edz-qxxdz; v.ri[10]=-qxzdx;
        v.ri[13]=v.ri[2]-dzqzz; v.ri[14]=v.ri[2]-dzqxx; v.ri[15]=-dxqxz;
        v.ri[16]=ee+eqzz+(v.ri[3]-ee)+qzzqzz; v.ri[17]=ee+eqzz+(v.ri[4]-ee)+qxxqzz;
        v.ri[18]=ee+eqxx+(v.ri[3]-ee)+qzzqxx; v.ri[19]=ee+eqxx+(v.ri[4]-ee)+qxxqxx;
        v.ri[20]=qxzqxz; v.ri[21]=ee+eqxx+(v.ri[4]-ee)+qxxqyy; v.ri[22]=.5*(qxxqxx-qxxqyy);
    }
finish:
    for (std::size_t i=2;i<23;++i) v.ri[i]*=ev_per_atomic_unit;
    v.core[2][1]=-static_cast<double>(b.core_charge)*v.ri[2]; v.core[3][1]=-static_cast<double>(b.core_charge)*v.ri[3]; v.core[4][1]=-static_cast<double>(b.core_charge)*v.ri[4];
    v.core[2][2]=-static_cast<double>(a.core_charge)*v.ri[5]; v.core[3][2]=-static_cast<double>(a.core_charge)*v.ri[11]; v.core[4][2]=-static_cast<double>(a.core_charge)*v.ri[12];
    return v;
}

struct Frame { std::array<double,3> x{},y{},z{}; std::array<std::array<double,8>,7> rot{}; };
Frame frame(const Atom& a, const Atom& b, const double r) {
    Frame f; f.x={{(a.x_angstrom-b.x_angstrom)/r,(a.y_angstrom-b.y_angstrom)/r,(a.z_angstrom-b.z_angstrom)/r}};
    if (std::abs(f.x[2])>.999999999999) { f.x[2]=std::copysign(1.0,f.x[2]); f.y={{0,1,0}}; f.z={{1,0,0}}; }
    else { const double z3=std::hypot(f.x[0],f.x[1]), q=1.0/z3; f.y={{-q*f.x[1]*std::copysign(1.0,f.x[0]),std::abs(q*f.x[0]),0}}; f.z={{-q*f.x[0]*f.x[2],-q*f.x[1]*f.x[2],z3}}; }
    int ij=0; for(int i=0;i<3;++i) for(int j=0;j<=i;++j){++ij; f.rot[static_cast<std::size_t>(ij)][1]=f.x[static_cast<std::size_t>(i)]*f.x[static_cast<std::size_t>(j)]; f.rot[static_cast<std::size_t>(ij)][2]=f.y[static_cast<std::size_t>(i)]*f.y[static_cast<std::size_t>(j)]+f.z[static_cast<std::size_t>(i)]*f.z[static_cast<std::size_t>(j)];}
    return f;
}

std::vector<double> rotate_w(ReppValues v, const Frame& initial, const bool heavy_a, const bool heavy_b) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::IntegralPreprocessing);
    auto f=initial; std::array<double,101>w{}; std::array<double,100>buf{};
    for(int i:{2,3,4,5,11,12})v.ri[static_cast<std::size_t>(i)]*=.5;
    w[1]=v.ri[1]*.25; if(!heavy_a&&!heavy_b)return {w[1]}; int kbuf=1;
    if(heavy_a){for(int i=1;i<=3;++i){w[static_cast<std::size_t>(i+1)]=v.core[2][1]*f.x[static_cast<std::size_t>(i-1)];buf[static_cast<std::size_t>(i)]=v.ri[2]*f.x[static_cast<std::size_t>(i-1)];}for(int k=1;k<=2;++k)for(int i=1;i<=6;++i){w[static_cast<std::size_t>(i+4)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*v.core[static_cast<std::size_t>(k+2)][1];buf[static_cast<std::size_t>(i+3)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*v.ri[static_cast<std::size_t>(k+2)];}for(int i:{4,6,9})buf[static_cast<std::size_t>(i)]*=.5;kbuf=10;}
    if(heavy_b){for(int i=2;i<=4;++i){w[static_cast<std::size_t>(i)]=v.core[2][2]*f.x[static_cast<std::size_t>(i-2)];buf[static_cast<std::size_t>(kbuf++)]=v.ri[5]*f.x[static_cast<std::size_t>(i-2)];}for(int k=1;k<=2;++k)for(int i=1;i<=6;++i){w[static_cast<std::size_t>(i+4)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*v.core[static_cast<std::size_t>(k+2)][2];buf[static_cast<std::size_t>(i+kbuf-1)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*v.ri[static_cast<std::size_t>(k+10)];}for(int i:{kbuf,kbuf+2,kbuf+5})buf[static_cast<std::size_t>(i)]*=.5;}
    if(!(heavy_a&&heavy_b)){for(std::size_t i=0;i<l1scat.size();++i)w[static_cast<std::size_t>(l1scat[i])]=buf[i+1];return {w.begin()+1,w.begin()+11};}
    for(int k=1;k<=2;++k) {
        for(int i=1;i<=6;++i)buf[static_cast<std::size_t>(i+18)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*v.ri[static_cast<std::size_t>(k+5)];
    }
    buf[25]=buf[20];
    buf[26]=buf[22];
    buf[27]=buf[23];
    int ij=0;for(int i=0;i<3;++i)for(int j=0;j<=i;++j){++ij;auto& rr=f.rot[static_cast<std::size_t>(ij)];rr[3]=f.x[static_cast<std::size_t>(i)]*f.y[static_cast<std::size_t>(j)]+f.x[static_cast<std::size_t>(j)]*f.y[static_cast<std::size_t>(i)];rr[4]=f.x[static_cast<std::size_t>(i)]*f.z[static_cast<std::size_t>(j)]+f.x[static_cast<std::size_t>(j)]*f.z[static_cast<std::size_t>(i)];rr[5]=f.y[static_cast<std::size_t>(i)]*f.z[static_cast<std::size_t>(j)]+f.y[static_cast<std::size_t>(j)]*f.z[static_cast<std::size_t>(i)];rr[6]=f.y[static_cast<std::size_t>(i)]*f.y[static_cast<std::size_t>(j)];rr[7]=f.z[static_cast<std::size_t>(i)]*f.z[static_cast<std::size_t>(j)];}
    for(int col=1;col<=7;++col)for(int row:{1,3,6})f.rot[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]*=.5;
    ij=2;for(int i=0;i<3;++i){w[static_cast<std::size_t>(ij)]=v.ri[8]*f.x[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+1)]=v.ri[9]*f.x[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+2)]=v.ri[10]*f.y[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+3)]=v.ri[10]*f.z[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+4)]=v.ri[13]*f.x[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+5)]=v.ri[14]*f.x[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+6)]=v.ri[15]*f.y[static_cast<std::size_t>(i)];w[static_cast<std::size_t>(ij+7)]=v.ri[15]*f.z[static_cast<std::size_t>(i)];ij+=8;}
    int ibuf=27,l=1;for(int j=1;j<=6;++j){for(int k=1;k<=4;++k){++l;for(int i=1;i<=6;++i)buf[static_cast<std::size_t>(i+ibuf)]+=f.rot[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*w[static_cast<std::size_t>(l)];}ibuf+=6;}
    ij=2;for(int i=1;i<=6;++i){w[static_cast<std::size_t>(ij)]=v.ri[16]*f.rot[static_cast<std::size_t>(i)][1]+v.ri[17]*f.rot[static_cast<std::size_t>(i)][2];w[static_cast<std::size_t>(ij+1)]=v.ri[18]*f.rot[static_cast<std::size_t>(i)][1];w[static_cast<std::size_t>(ij+2)]=v.ri[20]*f.rot[static_cast<std::size_t>(i)][3];w[static_cast<std::size_t>(ij+3)]=v.ri[20]*f.rot[static_cast<std::size_t>(i)][4];w[static_cast<std::size_t>(ij+4)]=v.ri[22]*f.rot[static_cast<std::size_t>(i)][5];w[static_cast<std::size_t>(ij+5)]=v.ri[19]*f.rot[static_cast<std::size_t>(i)][6]+v.ri[21]*f.rot[static_cast<std::size_t>(i)][7];w[static_cast<std::size_t>(ij+6)]=v.ri[19]*f.rot[static_cast<std::size_t>(i)][7]+v.ri[21]*f.rot[static_cast<std::size_t>(i)][6];ij+=7;}
    for(int col=0;col<6;++col)for(int row=0;row<6;++row){double sum=0;for(int k=0;k<7;++k)sum+=f.rot[static_cast<std::size_t>(row+1)][static_cast<std::size_t>(k+1)]*w[static_cast<std::size_t>(2+col*7+k)];buf[static_cast<std::size_t>(64+col*6+row)]=sum;}
    for(std::size_t i=0;i<l2scat.size();++i)w[static_cast<std::size_t>(l2scat[i])]=buf[i+1];
    return {w.begin()+1,w.end()};
}

void add_attraction(std::vector<double>& h, const std::size_t n, const std::size_t first,
                    const int orbitals, const int column, const ReppValues& v, const Frame& f) {
    h[first*n+first]+=v.core[1][static_cast<std::size_t>(column)]; if(orbitals==1)return;
    constexpr std::array<std::pair<int,int>,10> slots{{{0,0},{1,0},{1,1},{2,0},{2,1},{2,2},{3,0},{3,1},{3,2},{3,3}}};
    constexpr std::array<int,3> pslot{{1,3,6}}; for(int axis=0;axis<3;++axis){const double x=v.core[2][static_cast<std::size_t>(column)]*f.x[static_cast<std::size_t>(axis)];const auto [r,c]=slots[static_cast<std::size_t>(pslot[static_cast<std::size_t>(axis)])];h[(first+static_cast<std::size_t>(r))*n+first+static_cast<std::size_t>(c)]+=x;h[(first+static_cast<std::size_t>(c))*n+first+static_cast<std::size_t>(r)]+=x;}
    constexpr std::array<int,6> qslot{{2,4,5,7,8,9}};for(int k=0;k<6;++k){const double x=f.rot[static_cast<std::size_t>(k+1)][1]*v.core[3][static_cast<std::size_t>(column)]+f.rot[static_cast<std::size_t>(k+1)][2]*v.core[4][static_cast<std::size_t>(column)];const auto [r,c]=slots[static_cast<std::size_t>(qslot[static_cast<std::size_t>(k)])];h[(first+static_cast<std::size_t>(r))*n+first+static_cast<std::size_t>(c)]+=x;if(r!=c)h[(first+static_cast<std::size_t>(c))*n+first+static_cast<std::size_t>(r)]+=x;}
}

}  // namespace

ElectronicSystem build_electronic_system(std::span<const Atom> atoms, const CalculationOptions& options) {
    validate_cartesian_atoms(atoms);
    if(options.multiplicity!=1) {
        throw UnsupportedCalculationError("multiplicity " + std::to_string(options.multiplicity)
            + " is unsupported; only closed-shell singlet multiplicity 1 is supported");
    }
    ElectronicSystem s;
    std::int64_t electrons=-static_cast<std::int64_t>(options.molecular_charge);
    {
        AMSOLCPP_PERF_SCOPE(PerformanceStage::BasisAndParameters);
        s.atoms.assign(atoms.begin(),atoms.end());
        for(std::size_t i=0;i<atoms.size();++i){const auto& a=atoms[i];if(!is_supported_element(a.atomic_number))throw UnsupportedElementError(a.atomic_number,"AM1");if(!std::isfinite(a.x_angstrom)||!std::isfinite(a.y_angstrom)||!std::isfinite(a.z_angstrom))throw InputError("atom coordinates must be finite");const auto& p=element_parameters(a.atomic_number);s.parameters.push_back(&p);s.first_orbital.push_back(s.orbital_count);for(int k=0;k<p.basis_orbitals;++k)s.ao_to_atom.push_back(i);s.orbital_count+=static_cast<std::size_t>(p.basis_orbitals);electrons+=p.valence_electrons;}
    }
    if(electrons<=0||electrons%2!=0) {
        throw UnsupportedCalculationError("molecular charge "
            + std::to_string(options.molecular_charge) + " gives "
            + std::to_string(electrons)
            + " valence electrons; closed-shell AM1 requires a positive even count");
    }
    s.electron_count=static_cast<std::size_t>(electrons);
    s.occupied_orbitals=s.electron_count/2;
    if(s.occupied_orbitals>s.orbital_count) {
        throw UnsupportedCalculationError("electron count " + std::to_string(s.electron_count)
            + " for molecular charge " + std::to_string(options.molecular_charge)
            + " exceeds minimal-basis capacity " + std::to_string(2 * s.orbital_count));
    }
    const std::size_t n=s.orbital_count;
    {
        AMSOLCPP_PERF_SCOPE(PerformanceStage::OneElectronTerms);
        s.hcore.assign(n*n,0.0);for(std::size_t ai=0;ai<atoms.size();++ai){const auto& p=*s.parameters[ai];const auto first=s.first_orbital[ai];s.hcore[first*n+first]=p.uss;for(int k=1;k<p.basis_orbitals;++k)s.hcore[(first+static_cast<std::size_t>(k))*n+first+static_cast<std::size_t>(k)]=p.upp;}
    }
    for(std::size_t ai=1;ai<atoms.size();++ai)for(std::size_t aj=0;aj<ai;++aj){
        AMSOLCPP_PERF_INCREMENT(integral_pair_evaluations);
        const auto& pa=*s.parameters[ai];const auto& pb=*s.parameters[aj];const double r=distance(atoms[ai],atoms[aj]);if(r<1.0e-8)throw InputError("two atoms occupy the same Cartesian position");const auto ov=overlap_block(atoms[ai],atoms[aj],pa,pb);for(int u=0;u<pa.basis_orbitals;++u)for(int q=0;q<pb.basis_orbitals;++q){const double ba=u==0?pa.betas:pa.betap,bb=q==0?pb.betas:pb.betap;const std::size_t row=s.first_orbital[ai]+static_cast<std::size_t>(u),col=s.first_orbital[aj]+static_cast<std::size_t>(q);const double value=ov[static_cast<std::size_t>(u*pb.basis_orbitals+q)]*.5*(ba+bb);s.hcore[row*n+col]=value;s.hcore[col*n+row]=value;}
        // REPP uses the historical AMSOL bohr convention internally.
        auto rv=repp(r/bohr_radius_angstrom,pa,pb);const auto fr=frame(atoms[ai],atoms[aj],r);add_attraction(s.hcore,n,s.first_orbital[ai],pa.basis_orbitals,1,rv,fr);add_attraction(s.hcore,n,s.first_orbital[aj],pb.basis_orbitals,2,rv,fr);
        auto w=rotate_w(rv,fr,pa.basis_orbitals==4,pb.basis_orbitals==4);const std::size_t ni=static_cast<std::size_t>(pa.basis_orbitals*(pa.basis_orbitals+1)/2),nj=static_cast<std::size_t>(pb.basis_orbitals*(pb.basis_orbitals+1)/2);for(std::size_t row=0;row<ni;++row)if(product_pairs[row].first==product_pairs[row].second)for(std::size_t col=0;col<nj;++col)w[row*nj+col]*=2;for(std::size_t col=0;col<nj;++col)if(product_pairs[col].first==product_pairs[col].second)for(std::size_t row=0;row<ni;++row)w[row*nj+col]*=2;s.integral_blocks.push_back({ai,aj,static_cast<std::size_t>(pa.basis_orbitals),static_cast<std::size_t>(pb.basis_orbitals),std::move(w)});
        double scale=std::exp(-pa.alpha*r)+std::exp(-pb.alpha*r);if(pa.atomic_number+pb.atomic_number==8||pa.atomic_number+pb.atomic_number==9){if(pa.atomic_number==7||pa.atomic_number==8)scale+=(r-1)*std::exp(-pa.alpha*r);if(pb.atomic_number==7||pb.atomic_number==8)scale+=(r-1)*std::exp(-pb.alpha*r);}const double bare=static_cast<double>(pa.core_charge*pb.core_charge)*rv.ri[1];double scalar=0;for(std::size_t k=0;k<pa.gaussian_term_count;++k){const auto& g=pa.gaussian_terms[k];scalar+=g.amplitude*std::exp(-g.exponent*(r-g.displacement)*(r-g.displacement));}for(std::size_t k=0;k<pb.gaussian_term_count;++k){const auto& g=pb.gaussian_terms[k];scalar+=g.amplitude*std::exp(-g.exponent*(r-g.displacement)*(r-g.displacement));}s.nuclear_repulsion_ev+=bare+std::abs(scale*bare)+scalar*static_cast<double>(pa.core_charge*pb.core_charge)/r;
    }return s;
}

}  // namespace amsolcpp::detail
