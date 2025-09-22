#!/bin/bash
#
# Build script for combined stochastic sampling module
# Always work in msani_dev conda environment as per instructions
#

set -e  # Exit on any error

echo "=== Combined Stochastic Sampling Build Script ==="
echo ""

# Activate msani_dev conda environment
echo "Activating msani_dev conda environment..."
if ! conda info --envs | grep -q "msani_dev"; then
    echo "Error: msani_dev conda environment not found!"
    echo "Please create the environment first or check the name."
    exit 1
fi

eval "$(conda shell.bash hook)"
conda activate msani_dev

echo "Current conda environment: $CONDA_DEFAULT_ENV"
echo "Conda prefix: $CONDA_PREFIX"
echo ""

# Set up environment variables
export RDBASE="$CONDA_PREFIX"
export PYTHONPATH="$CONDA_PREFIX/lib/python3.11/site-packages:$PYTHONPATH"

# Clean previous build
echo "Cleaning previous build..."
rm -rf build/
mkdir -p build

# Configure with CMake
echo "Configuring with CMake..."
cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
    -DCMAKE_FIND_ROOT_PATH="$CONDA_PREFIX" \
    -DPython_EXECUTABLE="$CONDA_PREFIX/bin/python" \
    -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
    -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG" \
    -DBUILD_SHARED_LIBS=OFF

echo ""
echo "Building the module..."
make -j$(nproc) stochastic_sampling_combined

if [ $? -eq 0 ]; then
    echo ""
    echo "=== BUILD SUCCESSFUL ==="
    echo ""
    echo "Module location: $(pwd)/stochastic_sampling_combined$(python3-config --extension-suffix)"
    echo ""
    echo "Testing the module..."
    
    # Test import
    python -c "
try:
    import stochastic_sampling_combined
    print('✅ Module import successful!')
    print(f'Module version: {stochastic_sampling_combined.__version__}')
    print(f'Available functions: {[name for name in dir(stochastic_sampling_combined) if not name.startswith(\"_\")]}')
except Exception as e:
    print(f'❌ Module import failed: {e}')
    exit(1)
"
    cp *.so ..
    echo ""
    echo "=== MODULE READY FOR USE ==="
    echo ""
    echo "To use the module:"
    echo "1. Copy the .so file to your target machine (same OS)"
    echo "2. Import in Python: import stochastic_sampling_combined"
    echo "3. Use stochastic_sampling_combined.stochastic_sampling_discrete() or"
    echo "   stochastic_sampling_combined.stochastic_sampling_continuous()"
    echo ""
    echo "Module file: stochastic_sampling_combined$(python3-config --extension-suffix)"
    echo "Size: $(du -h stochastic_sampling_combined$(python3-config --extension-suffix) | cut -f1)"
    
else
    echo ""
    echo "=== BUILD FAILED ==="
    echo "Check the error messages above and ensure all dependencies are installed."
    exit 1
fi