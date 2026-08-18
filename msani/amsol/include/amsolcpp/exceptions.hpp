#pragma once

#include <stdexcept>
#include <string>

namespace amsolcpp {

class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InputError : public Error {
public:
    using Error::Error;
};

class UnsupportedElementError : public InputError {
public:
    UnsupportedElementError(int atomic_number, const std::string& component);
};

class UnsupportedCalculationError : public InputError {
public:
    using InputError::InputError;
};

class NumericalError : public Error {
public:
    using Error::Error;
};

class DiagonalizationError : public NumericalError {
public:
    using NumericalError::NumericalError;
};

class ScfConvergenceError : public NumericalError {
public:
    using NumericalError::NumericalError;
};

}  // namespace amsolcpp
