# For Windows

Download the precompiled executable AMSOL from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi) and renamed it as amsol7.1.exe. The program automatically detects if the OS is Windows and will use this for calculation of desolvation.


# For Linux 

1. Install notes

You need to download, compile & install AMSOL7.1 here. The executable should be named amsol7.1

The precompiled AMSOL should work fine if put it in this folder with the correct name.

2. Instruction on how to compile

Download the source code from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi)

Extract the zip file and in terminal:

```bash
cd amsol7.1
sed -i '' "s/set F77  = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -O -o'/set F77  = 'gfortran -c -finit-local-zero -fno-automatic -ffixed-line-length-72 -std=legacy -Iinclude -O -o'/g" amsol.compile && \
sed -i '' "s/set F77o = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -o'/set F77o = 'gfortran -c -finit-local-zero -fno-automatic -ffixed-line-length-72 -std=legacy -Iinclude -o'/g" amsol.compile && \
sed -i '' "s/set LD   = 'g77 -o'/set LD   = 'gfortran -ffixed-line-length-72 -o'/g" amsol.compile && \
sed -i 's/ OPEN(20,NAME=/ OPEN(20,FILE='/g new/amsol.f
sed -i 's/ OPEN(19,NAME=/ OPEN(19,FILE='/g new/amsol.f
csh amsol.compile<<EOF
man
linux
amsol7.1
s
EOF
chmod +x amsol7.1
```

Move amsol7.1 to this folder. Enjoy sanitizing molecules ;)

# For MacOS - Intel based

## Prerequisites

Install Homebrew if you don’t already have it:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Install the necessary build tools:

```bash
brew install gcc
brew install make
brew install tcsh
```

This provides gfortran, gcc, and the C-shell (csh) needed to run AMSOL’s legacy build script.

## Prepare the source code

Download the source code from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi)

Extract the zip file and in terminal:
```bash
cd amsol7.1 && \
sed -i '' "s/set F77  = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -O -o'/set F77  = 'gfortran -std=legacy -c -finit-local-zero -fno-automatic -Iinclude -O -o'/g" amsol.compile && \
sed -i '' "s/set F77o = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -o'/set F77o = 'gfortran -std=legacy -c -finit-local-zero -fno-automatic -Iinclude -o'/g" amsol.compile && \
sed -i '' "s/set LD   = 'g77 -o'/set LD   = 'gfortran -o'/g" amsol.compile && \
sed -i '' "s/ OPEN(20,NAME=/ OPEN(20,FILE=/" new/amsol.f && \
sed -i '' "s/ OPEN(19,NAME=/ OPEN(19,FILE=/" new/amsol.f && \
csh amsol.compile <<EOF
man
linux
amsol7.1_macos_x64
s
EOF
chmod +x amsol7.1_macos_x64
```

Then move `amsol7.1_macos_x64` to this folder. The program will detect automatically the architecture of the machine and run the software accordingly.

# For MacOS ARM64

## Prerequisites

Install Homebrew if you don’t already have it:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Install the necessary build tools:

```bash
brew install gcc
brew install make
brew install tcsh
```

## Prepare the source code

Download the source code from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi)

Extract the zip file then in terminal, use:

```bash
cd amsol7.1 && \
sed -i '' "s/set F77  = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -O -o'/set F77  = 'gfortran -std=legacy -c -finit-local-zero -fno-automatic -Iinclude -O -o'/g" amsol.compile && \
sed -i '' "s/set F77o = 'g77 -c -finit-local-zero -fno-automatic -Iinclude -o'/set F77o = 'gfortran -std=legacy -c -finit-local-zero -fno-automatic -Iinclude -o'/g" amsol.compile && \
sed -i '' "s/set LD   = 'g77 -o'/set LD   = 'gfortran -o'/g" amsol.compile && \
sed -i '' "s/ OPEN(20,NAME=/ OPEN(20,FILE=/" new/amsol.f && \
sed -i '' "s/ OPEN(19,NAME=/ OPEN(19,FILE=/" new/amsol.f && \
csh amsol.compile <<EOF
man
linux
amsol7.1_macos_arm64
s
EOF
chmod +x amsol7.1_macos_arm64
```

Then move `amsol7.1_macos_arm64` to this folder. The program will detect automatically the architecture of the machine and run the software accordingly.