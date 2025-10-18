// Boost throw_exception stub for Windows MSVC builds
// This provides the boost::throw_exception function that RDKit expects
// On Windows with MSVC, Boost headers declare but don't define this function

#include <exception>
#include <boost/throw_exception.hpp>

#ifdef _MSC_VER

namespace boost {

// Provide implementation of throw_exception for MSVC
#ifdef BOOST_NO_EXCEPTIONS
void throw_exception(std::exception const & e) {
    // If exceptions are disabled, terminate
    std::terminate();
}
#else
void throw_exception(std::exception const & e) {
    throw e;
}
#endif

} // namespace boost

#endif // _MSC_VER
