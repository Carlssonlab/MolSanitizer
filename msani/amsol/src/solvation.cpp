#include "science.hpp"
#include "performance_diagnostics.hpp"

#include <amsolcpp/exceptions.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace amsolcpp::detail {
namespace {

constexpr double pi = std::numbers::pi_v<double>;
constexpr double four_pi = 4.0 * pi;
constexpr double iofr = 1.4345;

struct SmParameter {
    int z; double surface_radius; double coulomb_radius; double water_sigma;
    double water_hrfn; double genorg_sigma_n; double genorg_hrfn_n;
};
constexpr std::array<SmParameter,11> sm_parameters{{
    {1,1.20,.91,99.48,0,39.40,0}, {6,1.70,1.78,113.37,-151.17,70.44,-112.19},
    {7,1.55,1.92,37.02,-270.26,42.07,-95.42}, {8,1.52,1.60,-148.98,-243.75,-24.56,-10.73},
    {9,1.47,1.50,44.97,0,1.57,0}, {14,2.10,2.10,68.60,0,-118.50,0},
    {15,1.80,2.40,-40.88,0,-29.90,0}, {16,1.80,2.15,-56.49,-34.36,-74.94,37.67},
    {17,1.75,2.13,-2.72,0,-36.15,0}, {35,1.85,2.31,-20.37,0,-47.41,0},
    {53,1.98,2.66,-18.82,0,-50.90,0}}};

const SmParameter& sm_parameter(const int z) {
    const auto found=std::find_if(sm_parameters.begin(),sm_parameters.end(),[z](const auto& p){return p.z==z;});
    if(found==sm_parameters.end())throw UnsupportedElementError(z,"SM5.42R");
    return *found;
}

std::vector<double> distances(std::span<const Atom> atoms) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::GeometryPreprocessing);
    const std::size_t n=atoms.size();std::vector<double>d(n*n,0);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j){const double value=std::hypot(atoms[i].x_angstrom-atoms[j].x_angstrom,atoms[i].y_angstrom-atoms[j].y_angstrom,atoms[i].z_angstrom-atoms[j].z_angstrom);if(value<=0||!std::isfinite(value))throw InputError("degenerate Cartesian geometry in SM5.42R");d[i*n+j]=d[j*n+i]=value;}
    return d;
}

std::pair<std::vector<double>,std::vector<double>> gauss_legendre(const int order) {
    std::vector<double>x(static_cast<std::size_t>(order)),w(static_cast<std::size_t>(order));const int half=(order+1)/2;
    for(int i=0;i<half;++i){double z=std::cos(pi*(static_cast<double>(i)+.75)/(static_cast<double>(order)+.5)),previous=0,derivative=0;do{previous=z;double p1=1,p2=0;for(int j=1;j<=order;++j){const double p3=p2;p2=p1;p1=((2.0*j-1)*z*p2-(j-1.0)*p3)/j;}derivative=order*(z*p1-p2)/(z*z-1);z=previous-p1/derivative;}while(std::abs(z-previous)>2e-15);x[static_cast<std::size_t>(i)]=-z;x[static_cast<std::size_t>(order-1-i)]=z;const double weight=2/((1-z*z)*derivative*derivative);w[static_cast<std::size_t>(i)]=w[static_cast<std::size_t>(order-1-i)]=weight;}return {x,w};
}

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

[[nodiscard]] double dot(const Vec3& left, const Vec3& right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

[[nodiscard]] Vec3 cross(const Vec3& left, const Vec3& right) {
    return {
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x
    };
}

[[nodiscard]] Vec3 scale(const Vec3& value, const double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

[[nodiscard]] Vec3 add(const Vec3& left, const Vec3& right) {
    return {left.x + right.x, left.y + right.y, left.z + right.z};
}

[[nodiscard]] double source_sqrt(const double value, const char* label) {
    if (value < 0.0) {
        if (value > -1.0e-12) {
            return 0.0;
        }
        throw NumericalError(std::string(label) + " became negative in DAREAL");
    }
    return std::sqrt(value);
}

[[nodiscard]] double source_acos_argument(const double value) {
    if (value < -1.0 || value > 1.0) {
        if (value < -1.0 - 1.0e-12 || value > 1.0 + 1.0e-12) {
            throw NumericalError("DAREAL acos argument outside source range");
        }
        return std::clamp(value, -1.0, 1.0);
    }
    return value;
}

[[nodiscard]] std::size_t matrix_index(
    const std::size_t row,
    const std::size_t column,
    const std::size_t dimension
) {
    return row * dimension + column;
}

[[nodiscard]] double polygon_angle(
    const std::size_t ia,
    const std::size_t ib,
    std::span<const double> ctheta,
    std::span<const double> inverse_sine,
    const std::size_t dimension
) {
    const auto index = [dimension](const std::size_t row, const std::size_t column) {
        return matrix_index(row, column, dimension);
    };
    const double phi1 = ctheta[index(ia, ib)]
        - ctheta[index(ia, ia)] * ctheta[index(ib, ib)];
    return std::acos(source_acos_argument(phi1 * inverse_sine[ia] * inverse_sine[ib]));
}

// Source-equivalent analytical spherical-cap union from AMSOL 7.1 DAREAL.
// Segment labels intentionally remain one-based to preserve the original
// traversal, connection ordering, and cusp-regularization decisions.
[[nodiscard]] double accessible_solid_angle(
    std::span<const Atom> atoms,
    std::span<const double> radii,
    std::span<const double> distances_by_pair,
    const std::size_t k
) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::SurfaceGeometry);
    constexpr double eps_i = 1.0e-11;
    constexpr double eps_c = 2.0e-3;
    constexpr double two_pi = 6.283185307179586476925286766558;
    constexpr double source_four_pi = two_pi + two_pi;
    constexpr double w0 = 3.0 * eps_c / 16.0;
    constexpr double w1 = -0.5;
    constexpr double w2 = 3.0 / (8.0 * eps_c);
    constexpr double w4 = -1.0 / (16.0 * eps_c * eps_c * eps_c);

    const std::size_t n = atoms.size();
    double rk = radii[k];
    if (rk <= 0.0) {
        return 0.0;
    }

    for (int restart_count = 0; restart_count < 8; ++restart_count) {
        const double eps_k = eps_i * rk;
        std::vector<std::size_t> crossing_atom(1, 0);
        std::vector<int> cusp_label(1, 0);
        for (std::size_t atom = 0; atom < n; ++atom) {
            if (atom == k || radii[atom] <= 0.0) {
                continue;
            }
            const double rik = distances_by_pair[atom * n + k];
            int label = -1;
            if (std::abs(rk + rik - radii[atom]) < eps_c) {
                label = 1;
            } else if (std::abs(rik + radii[atom] - rk) < eps_c) {
                label = 2;
            } else if (std::abs(rik - radii[atom] - rk) < eps_c) {
                label = 3;
            } else {
                if (rk + radii[atom] - rik < eps_k) {
                    continue;
                }
                if (rik - std::abs(rk - radii[atom]) < eps_k) {
                    if (rk <= radii[atom]) {
                        return 0.0;
                    }
                    continue;
                }
                label = 0;
            }
            crossing_atom.push_back(atom);
            cusp_label.push_back(label);
        }

        const std::size_t crossing_count = crossing_atom.size() - 1;
        if (crossing_count == 0) {
            return source_four_pi;
        }

        const std::size_t dimension = crossing_count + 1;
        const auto index = [dimension](const std::size_t row, const std::size_t column) {
            return matrix_index(row, column, dimension);
        };
        std::vector<Vec3> normals(dimension * dimension);
        std::vector<double> ctheta(dimension * dimension, 0.0);
        std::vector<double> stheta(dimension, 0.0);
        std::vector<unsigned char> connected(dimension * dimension, 0);

        const double rk_inverse = 0.5 / rk;
        const double rk_squared = rk * rk;
        for (std::size_t segment = 1; segment <= crossing_count; ++segment) {
            const std::size_t atom = crossing_atom[segment];
            const double rik = distances_by_pair[atom * n + k];
            normals[index(segment, segment)] = {
                (atoms[atom].x_angstrom - atoms[k].x_angstrom) / rik,
                (atoms[atom].y_angstrom - atoms[k].y_angstrom) / rik,
                (atoms[atom].z_angstrom - atoms[k].z_angstrom) / rik
            };
            const int label = cusp_label[segment];
            if (label == 0) {
                const double distance_inverse = 1.0 / rik;
                ctheta[index(segment, segment)] = rk_inverse
                    * (rik + (rk_squared - radii[atom] * radii[atom]) * distance_inverse);
                stheta[segment] = source_sqrt(
                    1.0 - ctheta[index(segment, segment)]
                        * ctheta[index(segment, segment)],
                    "DAREAL regular segment sine"
                );
            } else {
                double wx = 0.0;
                double cusp = 0.0;
                if (label == 1) {
                    wx = rik - radii[atom] + rk;
                    cusp = ((w4 * wx * wx + w2) * wx + w1) * wx + w0 + rik;
                } else if (label == 2) {
                    wx = rik + radii[atom] - rk;
                    cusp = ((w4 * wx * wx + w2) * wx + w1) * wx + w0 + rik;
                } else {
                    wx = rik - radii[atom] - rk;
                    cusp = ((-w4 * wx * wx - w2) * wx + w1) * wx - w0 + rik;
                }
                const double rik_inverse = 1.0 / cusp;
                double gamma = 0.0;
                if (label == 1) {
                    gamma = std::max(
                        1.0e-16,
                        ((rk + radii[atom]) * rik_inverse + 1.0)
                            * (cusp + rk - radii[atom]) * rik_inverse
                    );
                    ctheta[index(segment, segment)] = gamma - 1.0;
                } else {
                    gamma = std::max(
                        1.0e-16,
                        ((rk + radii[atom]) * rik_inverse - 1.0)
                            * (cusp - rk + radii[atom]) * rik_inverse
                    );
                    ctheta[index(segment, segment)] = 1.0 - gamma;
                }
                stheta[segment] = source_sqrt(
                    gamma * (2.0 - gamma),
                    "DAREAL cusp segment sine"
                );
            }
        }

        for (std::size_t i = 2; i <= crossing_count; ++i) {
            for (std::size_t j = 1; j < i; ++j) {
                if (connected[index(j, j)] != 0U) {
                    continue;
                }
                const double cisj = ctheta[index(i, i)] * stheta[j];
                const double sicj = stheta[i] * ctheta[index(j, j)];
                const double sisj = stheta[i] * stheta[j];
                ctheta[index(j, i)] = dot(normals[index(i, i)], normals[index(j, j)]);
                const double tij = ctheta[index(j, i)]
                    - ctheta[index(i, i)] * ctheta[index(j, j)];
                if (tij > sisj - eps_i * std::abs(cisj - sicj)) {
                    if (ctheta[index(j, j)] > ctheta[index(i, i)]) {
                        connected[index(j, j)] = 1U;
                    } else {
                        connected[index(i, i)] = 1U;
                        break;
                    }
                } else {
                    const double eps_ij = eps_i * (sicj + cisj);
                    bool pair_connected = false;
                    if (sicj + cisj >= 0.0) {
                        pair_connected = tij > eps_ij - sisj;
                    } else {
                        if (tij <= -sisj - eps_ij) {
                            return 0.0;
                        }
                        pair_connected = true;
                    }
                    connected[index(j, i)] = pair_connected ? 1U : 0U;
                    connected[index(i, j)] = pair_connected ? 1U : 0U;
                }
            }
        }

        double sliced_area = 0.0;
        std::vector<std::size_t> cluster_labels;
        std::vector<double> work_high(1, 0.0);
        std::vector<double> work_low(1, 0.0);
        for (std::size_t i = 1; i <= crossing_count; ++i) {
            if (connected[index(i, i)] != 0U) {
                continue;
            }
            bool isolated = true;
            for (std::size_t j = 1; j <= crossing_count; ++j) {
                if (connected[index(j, j)] == 0U && connected[index(j, i)] != 0U) {
                    isolated = false;
                    break;
                }
            }
            if (isolated) {
                sliced_area += 1.0 - ctheta[index(i, i)];
            } else {
                cluster_labels.push_back(i);
                work_high.push_back(ctheta[index(i, i)] + eps_i * stheta[i]);
                work_low.push_back(ctheta[index(i, i)] - eps_i * stheta[i]);
            }
        }
        sliced_area *= two_pi;
        if (cluster_labels.empty()) {
            return source_four_pi - sliced_area;
        }

        std::size_t free_intersections = 0;
        std::vector<std::vector<int>> connection_lists(dimension);
        bool restart = false;
        for (std::size_t order_i = 2; order_i <= cluster_labels.size(); ++order_i) {
            const std::size_t li = cluster_labels[order_i - 1];
            for (std::size_t order_j = 1; order_j < order_i; ++order_j) {
                const std::size_t lj = cluster_labels[order_j - 1];
                if (connected[index(lj, li)] == 0U) {
                    continue;
                }
                const double one_minus = 1.0
                    - ctheta[index(lj, li)] * ctheta[index(lj, li)];
                const double inverse_sin_squared = 1.0 / one_minus;
                const double aij = (ctheta[index(li, li)]
                    - ctheta[index(lj, lj)] * ctheta[index(lj, li)])
                    * inverse_sin_squared;
                const double bij = (ctheta[index(lj, lj)]
                    - ctheta[index(li, li)] * ctheta[index(lj, li)])
                    * inverse_sin_squared;
                const double cij = source_sqrt(
                    (1.0 - aij * ctheta[index(li, li)]
                        - bij * ctheta[index(lj, lj)]) * inverse_sin_squared,
                    "DAREAL free-intersection coefficient"
                );
                const Vec3 normal_cross = cross(
                    normals[index(li, li)], normals[index(lj, lj)]
                );
                const Vec3 base = add(
                    scale(normals[index(li, li)], aij),
                    scale(normals[index(lj, lj)], bij)
                );
                normals[index(li, lj)] = add(base, scale(normal_cross, cij));
                normals[index(lj, li)] = add(base, scale(normal_cross, -cij));

                bool free_ij = true;
                bool free_ji = true;
                for (std::size_t order_l = 1; order_l <= cluster_labels.size(); ++order_l) {
                    if (order_l == order_i || order_l == order_j) {
                        continue;
                    }
                    const std::size_t ll = cluster_labels[order_l - 1];
                    if (connected[index(ll, li)] == 0U
                        || connected[index(ll, lj)] == 0U) {
                        continue;
                    }
                    if (free_ji) {
                        const double check = dot(
                            normals[index(lj, li)], normals[index(ll, ll)]
                        );
                        if (check > work_high[order_l]) {
                            free_ji = false;
                        } else if (check >= work_low[order_l]) {
                            rk *= 1.0 + 4.0 * eps_i;
                            restart = true;
                            break;
                        }
                    }
                    if (free_ij) {
                        const double check = dot(
                            normals[index(li, lj)], normals[index(ll, ll)]
                        );
                        if (check > work_high[order_l]) {
                            free_ij = false;
                        } else if (check >= work_low[order_l]) {
                            rk *= 1.0 + 4.0 * eps_i;
                            restart = true;
                            break;
                        }
                    }
                    if (!free_ij && !free_ji) {
                        break;
                    }
                }
                if (restart) {
                    break;
                }
                if (!free_ij && !free_ji) {
                    continue;
                }
                ctheta[index(li, lj)] = ctheta[index(lj, li)];
                if (free_ji) {
                    ++free_intersections;
                    connection_lists[li].push_back(static_cast<int>(lj));
                    connection_lists[lj].push_back(-static_cast<int>(li));
                }
                if (free_ij) {
                    ++free_intersections;
                    connection_lists[li].push_back(-static_cast<int>(lj));
                    connection_lists[lj].push_back(static_cast<int>(li));
                }
            }
            if (restart) {
                break;
            }
        }
        if (restart) {
            continue;
        }
        if (free_intersections == 0) {
            return 0.0;
        }

        double polygon_area = 0.0;
        std::vector<double> inverse_sine(dimension, 0.0);
        for (const std::size_t li : cluster_labels) {
            auto& entries = connection_lists[li];
            const std::size_t nphi = entries.size();
            if (nphi == 0) {
                continue;
            }
            const int first = entries.front();
            const std::size_t lj = static_cast<std::size_t>(std::abs(first));
            const Vec3 connection_normal = first > 0
                ? normals[index(lj, li)]
                : normals[index(li, lj)];
            entries.front() = static_cast<int>(lj);
            inverse_sine[li] = 1.0 / stheta[li];
            const double c2i = ctheta[index(li, li)] * ctheta[index(li, li)];
            const Vec3 reference = cross(normals[index(li, li)], connection_normal);
            const bool polygon_on_left = dot(reference, normals[index(lj, lj)]) > 0.0;
            std::vector<std::pair<double, int>> angle_entries;
            angle_entries.reserve(nphi > 0 ? nphi - 1 : 0);
            for (std::size_t entry_index = 1; entry_index < nphi; ++entry_index) {
                const int entry = entries[entry_index];
                const std::size_t other = static_cast<std::size_t>(std::abs(entry));
                const Vec3 candidate = entry > 0
                    ? normals[index(other, li)]
                    : normals[index(li, other)];
                const double x_value = dot(reference, candidate);
                const double y_value = dot(candidate, connection_normal) - c2i;
                double angle = std::atan2(x_value, y_value);
                if (angle <= 0.0) {
                    angle += two_pi;
                }
                angle_entries.emplace_back(angle, static_cast<int>(other));
            }

            double odd_angle = 0.0;
            if (nphi == 2) {
                odd_angle = angle_entries.front().first;
                entries[1] = angle_entries.front().second;
            } else {
                std::sort(angle_entries.begin(), angle_entries.end());
                for (std::size_t entry_index = 0; entry_index < angle_entries.size(); ++entry_index) {
                    entries[entry_index + 1] = angle_entries[entry_index].second;
                }
                odd_angle = angle_entries.front().first;
                for (std::size_t fortran_j = 3; fortran_j < nphi; fortran_j += 2) {
                    odd_angle += angle_entries[fortran_j - 1].first
                        - angle_entries[fortran_j - 2].first;
                }
            }
            const double even_angle = two_pi - odd_angle;
            const double cap = 1.0 - ctheta[index(li, li)];
            if (polygon_on_left) {
                polygon_area += odd_angle;
                sliced_area += even_angle * cap;
            } else {
                polygon_area += even_angle;
                sliced_area += odd_angle * cap;
                std::rotate(entries.begin(), entries.begin() + 1, entries.end());
            }
        }

        std::size_t polygon_count = 0;
        for (const std::size_t li : cluster_labels) {
            auto& entries = connection_lists[li];
            for (std::size_t entry_index = 1; entry_index < entries.size(); entry_index += 2) {
                if (entries[entry_index] == 0) {
                    continue;
                }
                std::size_t ia = static_cast<std::size_t>(entries[entry_index - 1]);
                std::size_t ib = li;
                entries[entry_index] = 0;
                double phi = polygon_angle(ia, ib, ctheta, inverse_sine, dimension);
                polygon_area += phi;
                while (true) {
                    bool found = false;
                    auto& ia_entries = connection_lists[ia];
                    for (std::size_t list_index = 1;
                         list_index < ia_entries.size(); list_index += 2) {
                        if (ia_entries[list_index] != static_cast<int>(ib)) {
                            continue;
                        }
                        const std::size_t old_ib = ib;
                        ia_entries[list_index] = 0;
                        ib = ia;
                        ia = static_cast<std::size_t>(ia_entries[list_index - 1]);
                        if (ia != old_ib) {
                            phi = polygon_angle(
                                ia, ib, ctheta, inverse_sine, dimension
                            );
                        } else {
                            ib = li;
                            ia = static_cast<std::size_t>(entries[entry_index - 1]);
                        }
                        polygon_area += phi;
                        found = true;
                        break;
                    }
                    if (!found) {
                        break;
                    }
                }
                ++polygon_count;
            }
        }

        polygon_area += (static_cast<double>(polygon_count)
            - static_cast<double>(free_intersections)) * two_pi;
        return source_four_pi - sliced_area - std::fmod(polygon_area, source_four_pi);
    }
    throw NumericalError("DAREAL regularization restart did not converge");
}

int element_type(const int z) {
    switch(z){case 1:return 1;case 6:return 2;case 7:return 3;case 8:return 4;case 9:return 5;case 16:return 6;case 17:return 7;case 35:return 8;case 15:return 9;case 53:return 10;case 14:return 11;default:return 0;}
}
constexpr std::array<std::array<double,11>,11> rkk{{
{{0,1.55,1.55,1.55,0,2.14,0,0,0,0,0}},{{1.55,1.84,1.84,1.84,1.84,2.20,2.10,2.30,2.20,2.60,0}},
{{1.55,1.84,1.85,1.50,0,0,0,0,0,0,0}},{{1.55,1.84,1.50,2.75,0,0,0,0,2.10,0,2.10}},
{{0,1.84,0,0,0,0,0,0,0,0,0}},{{2.14,2.20,0,0,0,2.75,0,0,2.50,0,0}},
{{0,2.10,0,0,0,0,0,0,0,0,0}},{{0,2.30,0,0,0,0,0,0,0,0,0}},
{{0,2.20,0,2.10,0,2.50,0,0,0,0,0}},{{0,2.60,0,0,0,0,0,0,0,0,0}},
{{0,0,0,2.10,0,0,0,0,0,0,0}}}};
double rkk2(const int a,const int b){if(a==2&&b==2)return 1.27;if((a==2&&b==3)||(a==3&&b==2))return 1.225;if(a==2&&b==4)return 1.33;if(a==4&&b==1)return 1.55;if(a==4&&b==2)return 1.33;return 0;}
double tkk(const double d,const double cutoff,const double delta,const double factor=1){return cutoff<=0?0:std::exp(delta/(d-cutoff))*factor;}

struct Special {double cc4,cc5,oc,oo,nc,on,ss,cn2,nc2,hnn,hoh,op,sp,ntriple;};
Special special_parameters(const Solvent solvent){if(solvent==Solvent::Water)return{-67.68,16.10,231.29,41.82,-69.62,240.69,49.21,35.02,-238.94,-253.16,534.56,149.34,503.51,-3.55};return{-72.83*iofr,5.41*iofr,45.06*iofr,-30.49*iofr,-15.34*iofr,36.18*iofr,16.77*iofr,-69.86*iofr,0,-207.45*iofr,226.66*iofr,54.90*iofr,282.36*iofr,-14.41*iofr};}

double special_sigma(const std::size_t atom,std::span<const Atom> atoms,std::span<const double>d,
                     const Solvent solvent) {
    const std::size_t n=atoms.size();const int z=atoms[atom].atomic_number;const auto sp=special_parameters(solvent);const auto dist=[&](std::size_t i,std::size_t j){return d[i*n+j];};const auto hrfn=[solvent](int zz){const auto&p=sm_parameter(zz);return solvent==Solvent::Water?p.water_hrfn:p.genorg_hrfn_n*iofr;};double value=0;
    if(z==1){const double delta=.30;for(std::size_t j=0;j<n;++j){const int zj=atoms[j].atomic_number;if(zj!=6&&zj!=7&&zj!=8&&zj!=16)continue;const int tj=element_type(zj);const double cutoff=rkk[0][static_cast<std::size_t>(tj-1)]+delta;if(dist(j,atom)>=cutoff)continue;const double f0=tkk(dist(j,atom),cutoff,delta);value+=f0*hrfn(zj);if(zj==7&&sp.hnn!=0&&f0>0){const double nested_cutoff=rkk[static_cast<std::size_t>(tj-1)][2]+delta;double nested=0;for(std::size_t k=0;k<n;++k)if(j!=k&&atoms[k].atomic_number==7&&dist(k,j)<nested_cutoff)nested+=tkk(dist(k,j),nested_cutoff,delta,sp.hnn);value+=f0*nested;}else if(zj==8&&sp.hoh!=0&&f0>0){const double nested_cutoff=rkk[static_cast<std::size_t>(tj-1)][0]+delta;double nested=0;for(std::size_t k=0;k<n;++k)if(atom!=k&&atoms[k].atomic_number==1&&dist(k,j)<nested_cutoff)nested+=tkk(dist(k,j),nested_cutoff,delta,sp.hoh);value+=f0*nested;}}return value;}
    if(z==6){double delta=.30,cutoff=rkk[1][1]+delta;for(std::size_t j=0;j<n;++j)if(j!=atom&&atoms[j].atomic_number==6&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.cc4);delta=.07;cutoff=rkk2(2,2)+delta;for(std::size_t j=0;j<n;++j)if(j!=atom&&atoms[j].atomic_number==6&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.cc5);delta=.30;cutoff=rkk[1][2]+delta;double sum=0;for(std::size_t j=0;j<n;++j)if(atoms[j].atomic_number==7&&dist(j,atom)<cutoff)sum+=tkk(dist(j,atom),cutoff,delta);return value+sp.cn2*sum*sum;}
    if(z==7){const int carbon=2;double delta=.30,cutoff=rkk[2][1]+delta,sum1=0,sum2=0;for(std::size_t j=0;j<n;++j)if(atoms[j].atomic_number==6&&dist(j,atom)<cutoff){const double q=tkk(dist(j,atom),cutoff,delta);double neighbors=0,oxygen=0;for(std::size_t k=0;k<n;++k){if(k==atom||k==j)continue;const int tk=element_type(atoms[k].atomic_number);const double nested_cutoff=rkk[static_cast<std::size_t>(tk-1)][static_cast<std::size_t>(carbon-1)]+delta;if(nested_cutoff>delta&&dist(k,j)<nested_cutoff){const double x=tkk(dist(k,j),nested_cutoff,delta);neighbors+=x;if(atoms[k].atomic_number==8)oxygen+=x;}}sum1+=q*neighbors*neighbors;sum2+=q*oxygen;}value=std::pow(sum1,1.3)*sp.nc+sum2*sp.nc2;delta=.065;cutoff=rkk2(3,2)+delta;for(std::size_t j=0;j<n;++j)if(j!=atom&&atoms[j].atomic_number==6&&cutoff>delta&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.ntriple);return value;}
    if(z==8){double delta=.10,cutoff=rkk2(4,2)+delta;for(std::size_t j=0;j<n;++j)if(atoms[j].atomic_number==6&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.oc);delta=.30;cutoff=rkk[3][2]+delta;for(std::size_t j=0;j<n;++j)if(atoms[j].atomic_number==7&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.on);cutoff=rkk[3][3]+delta;double sum=0;int count=0;for(std::size_t j=0;j<n;++j)if(j!=atom&&atoms[j].atomic_number==8&&dist(j,atom)<cutoff){sum+=tkk(dist(j,atom),cutoff,delta);++count;}if(count&&sum>0)value+=sp.oo*std::exp(-1/(sum/.4));cutoff=rkk[3][8]+delta;for(std::size_t j=0;j<n;++j)if(atoms[j].atomic_number==15&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,sp.op);return value;}
    if(z==16){const double delta=.30;for(const auto& pair:std::array<std::pair<int,double>,2>{{{15,sp.sp},{16,sp.ss}}}){const int other_type=element_type(pair.first);const double cutoff=rkk[5][static_cast<std::size_t>(other_type-1)]+delta;for(std::size_t j=0;j<n;++j)if(j!=atom&&atoms[j].atomic_number==pair.first&&dist(j,atom)<cutoff)value+=tkk(dist(j,atom),cutoff,delta,pair.second);}return value;}
    return 0;
}

}  // namespace

SolvationModel build_solvation_model(std::span<const Atom> atoms,const Solvent solvent) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::SolvationOverall);
    if(atoms.empty())throw InputError("SM5.42R requires at least one atom");
    const std::size_t n=atoms.size();const auto d=distances(atoms);SolvationModel model;model.solvent=solvent;model.dielectric=solvent==Solvent::Water?78.3:2.06;model.coulomb_radii.reserve(n);model.surface_radii.reserve(n);model.base_sigma.reserve(n);
    for(const auto&a:atoms){const auto&p=sm_parameter(a.atomic_number);model.coulomb_radii.push_back(p.coulomb_radius);model.surface_radii.push_back(p.surface_radius);model.base_sigma.push_back(solvent==Solvent::Water?p.water_sigma:p.genorg_sigma_n*iofr);}
    const int nt=std::clamp(static_cast<int>(10.0+.2*static_cast<double>(n)),4,16);const auto radial=gauss_legendre(nt);model.born_radii.resize(n);
    for(std::size_t k=0;k<n;++k){auto radii=model.coulomb_radii;const double rinf=radii[k];double target=0,sum36=0;for(std::size_t i=0;i<n;++i){const double v=d[i*n+k]+radii[i];target=std::max(target,v);if(i!=k)sum36+=std::pow(v,36);}if(target/rinf<1.00001){model.born_radii[k]=rinf;continue;}const double rsupln=std::log(sum36)/36.0,rinfln=std::log(rinf),ascal=.5*(rsupln-rinfln),bscal=.5*(rsupln+rinfln),rsup=std::exp(rsupln);double integral=four_pi/rsup;for(std::size_t q=0;q<radial.first.size();++q){const double x=std::exp(radial.first[q]*ascal+bscal);radii[k]=x;integral+=accessible_solid_angle(atoms,radii,d,k)*radial.second[q]*ascal/x;}model.born_radii[k]=four_pi/integral;}
    model.fgb.assign(n*n,0);const double coefficient=14.3966*(1.0-1.0/model.dielectric);for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j){const bool hc=(atoms[i].atomic_number==1&&atoms[j].atomic_number==6)||(atoms[i].atomic_number==6&&atoms[j].atomic_number==1);const double dkk=hc?4.2:3.9,aij=model.born_radii[i]*model.born_radii[j],rij=d[i*n+j];const double value=coefficient/std::sqrt(rij*rij+aij*std::exp(-rij*rij/(dkk*aij)));model.fgb[i*n+j]=model.fgb[j*n+i]=value;}
    model.areas.resize(n);for(std::size_t k=0;k<n;++k)model.areas[k]=accessible_solid_angle(atoms,model.surface_radii,d,k)*model.surface_radii[k]*model.surface_radii[k];complete_cds(model,atoms);return model;
}

void complete_cds(SolvationModel& model,std::span<const Atom> atoms) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::Cds);
    const auto d=distances(atoms);model.special_sigma.resize(atoms.size());for(std::size_t i=0;i<atoms.size();++i)model.special_sigma[i]=special_sigma(i,atoms,d,model.solvent);if(model.solvent==Solvent::Hexadecane){double area=0;for(double x:model.areas)area+=x;model.large_surface_contribution=area*(.3871*38.93)/1000.0;}
}

}  // namespace amsolcpp::detail
