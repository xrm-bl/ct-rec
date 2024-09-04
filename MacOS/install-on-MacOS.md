# Install on MacOS
- You should make some modifications on Makefile and scripts for MacOS (both M1 and intel chips).
- Some essential points for modifications are described here.
- You can also refer modified src and bin files in this directory.
- This instruction and files were made by Macbook Air (2024, Apple M3 chip (8 cores CPU, 10 cores GPU, 16 cores Neural Engine), 24GB unified memory, OS: Sonoma 14.5).

## Step 0: Check your architecture
- Launch Terminal.
- `uname -m`
  - `arm64`: You're using M1 chip
  - `x86_64`: You're using Intel chip (or emulating intel on M1 chip)
- For M1 chip, you can switch `arm64` and `x86_64` by Rosetta2
  - Right-click on Terminal icon, click "Get summary info"("情報を見る").
  - Check "Open using Rosetta"("Rosettaを使用して開く") to launch with intel.
  - Uncheck "Open using Rosetta"("Rosettaを使用して開く") to launch with M1.
  - Don't forget to reboot launched Terminal after switching.

## Step 1: Prerequisites for MacOS
### Brew
- It's recommended to use brew for installing packages. https://brew.sh/ja/
### gcc
- Default compiler (clang) will cause errors. Instead, we use gcc.
- Install gcc ver.13.
```
brew install gcc@13
```
- Make symbolic links for gcc@13.   

For M1(`arm64`),
```
ln -s /opt/homebrew/bin/gcc-13 /usr/local/bin/gcc
```
For intel(`x86_64`),   
```
ln -s /usr/local/bin/gcc-13 /usr/local/bin/gcc
```
- Reboot Terminal.
### Libtiff
- Install libtiff.
```
brew install libtiff
```
- This instruction was made with libtiff ver. 4.6.0.

## Step 2: Modify scripts
### ct-rec/src/MakefileCPU
- Edit the first part of the file.   

For M1(`arm64`),
```.mk
INCLUDE	=-I/opt/homebrew/include
LIB	=-L/opt/homebrew/lib
BIN	=../bin

CC	=gcc -O3 -D_GNU_SOURCE $(INCLUDE) $(LIB)
CC0	=gcc -O0 -D_GNU_SOURCE $(INCLUDE) $(LIB)
```
For intel(`x86_64`), 
```.mk
INCLUDE	=-I/usr/local/include
LIB	=-L/usr/local/lib
BIN	=../bin

CC	=gcc -O3 -D_GNU_SOURCE $(INCLUDE) $(LIB)
CC0	=gcc -O0 -D_GNU_SOURCE $(INCLUDE) $(LIB)
```
- Edit `ofct_DO` part. (`ofct_DO.c` -> `oct_DO.c`)
```.mk
ofct_DO:	error.c rhp.h rhp.c msd.h msd.c sif_f.h $(SIF_F) oct_DO.c sif_f.h
	$(CC) oct_DO.c rhp.c msd.c $(ESF) -DONLY_CT_VIEWS $(LP) -lm -o ofct_DO
```
### ct-rec/src/oct_DO.c
- Add declaration in function `Scan`. (`HiPic		hp;`)
```oct_DO.c
static void	Scan(int Ox1,int Oy1,int Ox2,int Oy2,FOM *D,FOM *S,int m)
{
	int	y,x,
		X=Ox1,
		Y=Oy1;
	double	d=(*D);
	HiPic		hp;  <----Add like this
```

### ct-rec/src/rec_stk.c
- Comment out (or delete) `#include <malloc.h>`.
```rec_stk.c
//	program cthist
//	usage 'bmake > output_file'

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// #include <malloc.h>    <----This line
```

## Step 3: Make and install
- Make and install with `MakefileCPU`.
```.sh
cd ct-rec/src
make -f MakefileCPU all
mkdir ../bin
make -f MakefileCPU install
sudo cp ../bin/* /usr/local/bin
```
