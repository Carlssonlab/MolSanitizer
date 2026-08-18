#include <amsolcpp/elements.hpp>

#include <amsolcpp/exceptions.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <sstream>
#include <string>

namespace amsolcpp {
namespace {

using G = GaussianCoreTerm;

// AM1 parameters transcribed independently from AMSOL 7.1 PARAM.i and checked
// against the active pyAMSOL audit tables. Energies are eV and distances use
// the original AMSOL AM1 parameter conventions.
constexpr std::array<ElementParameters, 11> parameters{{
    {1, "H", 1, 1, 1,
     -11.3964270, 0.0, -6.1737870, 0.0, 1.1880780, 0.0, 2.8823240,
     -11.3964270, 12.8480000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.4721793, 0.4721793, 0.4721793,
     {G{0.1227960, 5.0000000, 1.2000000}, G{0.0050900, 5.0000000, 1.8000000},
      G{-0.0183360, 2.0000000, 2.1000000}, G{}}, 3},
    {6, "C", 4, 4, 4,
     -52.0286580, -39.6142390, -15.7157830, -7.7192830, 1.8086650, 1.6851160, 2.6482740,
     -120.8157940, 12.2300000, 11.4700000, 11.0800000, 9.8400000, 2.4300000,
     0.8236736, 0.7268015, 0.4494671, 0.6082946, 0.6423492,
     {G{0.0113550, 5.0000000, 1.6000000}, G{0.0459240, 5.0000000, 1.8500000},
      G{-0.0200610, 5.0000000, 2.0500000}, G{-0.0012600, 5.0000000, 2.6500000}}, 4},
    {7, "N", 5, 5, 4,
     -71.8600000, -57.1675810, -20.2991100, -18.2386660, 2.3154100, 2.1579400, 2.9472860,
     -202.4077430, 13.5900000, 12.6600000, 12.9800000, 11.5900000, 3.1400000,
     0.6433247, 0.5675528, 0.4994487, 0.7820840, 0.7883498,
     {G{0.0252510, 5.0000000, 1.5000000}, G{0.0289530, 5.0000000, 2.1000000},
      G{-0.0058060, 2.0000000, 2.4000000}, G{}}, 3},
    {8, "O", 6, 6, 4,
     -97.8300000, -78.2623800, -29.2727730, -29.2727730, 3.1080320, 2.5240390, 4.4553710,
     -316.0995200, 15.4200000, 14.4800000, 14.5200000, 12.9800000, 3.9400000,
     0.4988896, 0.4852322, 0.5667034, 0.9961066, 0.9065223,
     {G{0.2809620, 5.0000000, 0.8479180}, G{0.0814300, 7.0000000, 1.4450710}, G{}, G{}}, 2},
    {9, "F", 7, 7, 4,
     -136.1055790, -104.8898850, -69.5902770, -27.9223600, 3.7700820, 2.4946700, 5.5178000,
     -482.2905830, 16.9200000, 17.2500000, 16.7100000, 14.9100000, 4.8300000,
     0.4145203, 0.4909446, 0.6218302, 1.2088792, 0.9449355,
     {G{0.2420790, 4.8000000, 0.9300000}, G{0.0036070, 4.6000000, 1.6600000}, G{}, G{}}, 2},
    {14, "Si", 4, 4, 4,
     -33.9536220, -28.9347490, -3.7848520, -1.9681230, 1.8306970, 1.2849530, 2.2578160,
     -79.0017420, 9.8200000, 8.3600000, 7.3100000, 6.5400000, 1.3200000,
     1.1631107, 1.3022422, 0.3608967, 0.3829813, 0.3712106,
     {G{0.2500000, 9.0000000, 0.9114530}, G{0.0615130, 5.0000000, 1.9955690},
      G{0.0207890, 5.0000000, 2.9906100}, G{}}, 3},
    {15, "P", 5, 5, 4,
     -42.0298630, -34.0307090, -6.3537640, -6.5907090, 1.9812800, 1.8751500, 2.4553220,
     -124.4368355, 11.5600050, 5.2374490, 7.8775890, 7.3076480, 0.7792380,
     1.0452022, 0.8923660, 0.4248440, 0.3275319, 0.4386854,
     {G{-0.0318270, 6.0000000, 1.4743230}, G{0.0184700, 7.0000000, 1.7793540},
      G{0.0332900, 9.0000000, 3.0065760}, G{}}, 3},
    {16, "S", 6, 6, 4,
     -56.6940560, -48.7170490, -3.9205660, -7.9052780, 2.3665150, 1.6672630, 2.4616480,
     -191.7321930, 11.7863290, 8.6631270, 10.0393080, 7.7816880, 2.5321370,
     0.9004265, 1.0036329, 0.4331617, 0.5907115, 0.6454943,
     {G{-0.5091950, 4.5936910, 0.7706650}, G{-0.0118630, 5.8657310, 1.5033130},
      G{0.0123340, 13.5573360, 2.0091730}, G{}}, 3},
    {17, "Cl", 7, 7, 4,
     -111.6139480, -76.6401070, -24.5946700, -14.6372160, 3.6313760, 2.0767990, 2.9193680,
     -372.1984310, 15.0300000, 13.1600000, 11.3000000, 9.9700000, 2.4200000,
     0.5406286, 0.8057208, 0.5523705, 0.7693200, 0.6133369,
     {G{0.0942430, 4.0000000, 1.3000000}, G{0.0271680, 4.0000000, 2.1000000}, G{}, G{}}, 2},
    {35, "Br", 7, 7, 4,
     -104.6560630, -74.9300520, -19.3998800, -8.9571950, 3.0641330, 2.0383330, 2.5765460,
     -352.3142087, 15.0364395, 13.0346824, 11.2763254, 9.8544255, 2.4558683,
     0.8458104, 1.0407133, 0.5526071, 0.6024598, 0.5307555,
     {G{0.0666850, 4.0000000, 1.5000000}, G{0.0255680, 4.0000000, 2.3000000}, G{}, G{}}, 2},
    {53, "I", 7, 7, 4,
     -103.5896630, -74.4299970, -8.4433270, -6.3234050, 2.1028580, 2.1611530, 2.2994240,
     -346.8642857, 15.0404486, 13.0565580, 11.1477837, 9.9140907, 2.4563820,
     1.4878778, 1.1887388, 0.5527544, 0.4497523, 0.4631775,
     {G{0.0043610, 2.3000000, 1.8000000}, G{0.0157060, 3.0000000, 2.2400000}, G{}, G{}}, 2},
}};

constexpr std::array<int, parameters.size()> atomic_numbers{{1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}};

std::string canonical_symbol(std::string_view symbol) {
    std::string result;
    result.reserve(symbol.size());
    for (const char c : symbol) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
            result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        }
    }
    if (!result.empty()) {
        result.front() = static_cast<char>(std::toupper(static_cast<unsigned char>(result.front())));
    }
    return result;
}

}  // namespace

UnsupportedElementError::UnsupportedElementError(
    const int atomic_number_value,
    const std::string& component
) : InputError([&] {
    std::ostringstream message;
    message << "element with atomic number " << atomic_number_value
            << " is unsupported by " << component
            << "; supported atomic numbers are 1,6,7,8,9,14,15,16,17,35,53";
    return message.str();
}()) {}

const ElementParameters& element_parameters(const int atomic_number_value) {
    const auto found = std::find_if(parameters.begin(), parameters.end(), [&](const auto& item) {
        return item.atomic_number == atomic_number_value;
    });
    if (found == parameters.end()) {
        throw UnsupportedElementError(atomic_number_value, "AM1 parameter lookup");
    }
    return *found;
}

const ElementParameters& element_parameters(const std::string_view symbol) {
    const auto normalized = canonical_symbol(symbol);
    const auto found = std::find_if(parameters.begin(), parameters.end(), [&](const auto& item) {
        return item.symbol == normalized;
    });
    if (found == parameters.end()) {
        throw InputError("unsupported element symbol '" + std::string(symbol) + "'");
    }
    return *found;
}

bool is_supported_element(const int atomic_number_value) noexcept {
    return std::find(atomic_numbers.begin(), atomic_numbers.end(), atomic_number_value) != atomic_numbers.end();
}

int atomic_number(const std::string_view symbol) {
    return element_parameters(symbol).atomic_number;
}

std::string_view element_symbol(const int atomic_number_value) {
    return element_parameters(atomic_number_value).symbol;
}

std::span<const int> supported_atomic_numbers() noexcept {
    return atomic_numbers;
}

}  // namespace amsolcpp
