# For Windows

Download the precompiled executable AMSOL from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi) and renamed it as amsol7.1.exe. The program automatically detects if the OS is Windows and will use this for calculation of desolvation.


# For Linux 

1. Install notes

You need to download, compile & install AMSOL7.1 here. The executable should be named amsol7.1

The precompiled AMSOL should work fine if put it in this folder with the correct name.

2. Instruction on how to compile

Download the source code from [here](https://comp.chem.umn.edu/sds/amsol/amsol.cgi)

Extract the zip file and go to the AMSOL7.1 folder.

In lines 320-323 of amsol.compile, use your editor of interest, replace to:

```bash
set F77  = 'gfortran -c -finit-local-zero -fno-automatic -ffixed-line-length-72 -std=legacy -Iinclude -O -o'
set F77o = 'gfortran -c -finit-local-zero -fno-automatic -ffixed-line-length-72 -std=legacy -Iinclude -o'
set LD   = 'gfortran -ffixed-line-length-72 -o'
```

Save and return to the folder.
Run two modifications in command line:

```bash
sed -i 's/ OPEN(20,NAME=/ OPEN(20,FILE='/g new/amsol.f
sed -i 's/ OPEN(19,NAME=/ OPEN(19,FILE='/g new/amsol.f
```

Now it is able to compile AMSOL using (please read through the right below option before starting):
```bash
csh amsol.compile
```

Choose: man -> linux -> amsol7.1 -> s

Now the amsol is available. Don't forget to add the executable permission to such file:

```bash
chmod +x amsol7.1
```

Move amsol7.1 to this folder. Enjoy sanitizing molecules ;)
