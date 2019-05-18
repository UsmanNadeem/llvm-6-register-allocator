Register Allocator
=============================

This is a register allocator built (for a course project) on LLVM 6 that does not use the LLVM's existing RegAlloc interface, instead it is built using classes like TargetRegisterInfo, TargetRegisterClass, MachineFrameInfo etc.

It implements an liveness analysis based on the SSA form for virtual register. For physical registers the analysis is done at a basic block level beacause they are not live outside a BB. Spilling on the stack is done without splitting the live ranges, i.e. for the complete live range. The spilling cost/benefit ratio is calculated using the number of uses+defs and the number of edges in the interference graph. 

At this point, only small compiled programs work properly. Large programs either crash (during exection) or give incorrect output. There is probably some problem in the liveness analysis but I dont have time to debug and fix it. 


### Step I: Installing dependencies

1. `sudo apt-get install build-essential curl libcap-dev git cmake libncurses5-dev python-minimal python-pip unzip autoconf binutils-gold binutils-dev gcc-multilib zlib1g-dev`

2. `CMake 3.4.3` or higher:
	a. `wget https://cmake.org/files/v3.5/cmake-3.5.2.tar.gz`
	
    b. `tar xf cmake-3.5.2.tar.gz`
	
    c. `cd cmake-3.5.2`
	
    d. `cmake .`
	
    e. `make`

    f. `sudo make install`

3. Additionally some more LLVM dependencies might be required: `https://llvm.org/docs/GettingStarted.html#requirements`

### Step II: Clone and Build LLVM+allocator
1. Go to the LLVM folder downloaded from github or clone from this repo: `git clone https://github.com/UsmanNadeem/llvm-6-register-allocator.git`

2. `mkdir llvmbuild; cd llvmbuild`  --> make sure that the LLVM source folder is at the path `<../llvm-6-register-allocator>`
		
3. `cmake -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_RTTI=OFF -DLLVM_TARGETS_TO_BUILD="X86" -DLLVM_INCLUDE_TESTS=OFF  -DLLVM_BINUTILS_INCDIR=~/binutils/include  ../llvm-6-register-allocator`

4. `make llc clang clang++`

### Step III: Run the allocator benchmark
1. Compile to LLVM IR: `llvmbuild/bin/clang++ test.c -emit-llvm -S -o test.ll`

2. Use llc to create a text assembly file: `llvmbuild/bin/llc --custom-alloc=1 test.ll -o test.asm`. The flag to enable this allocator is `--custom-alloc=1`

3. Use clang to compile the assembly file to an executable: `llvmbuild/bin/clang++ test.asm -o test`

4. Run the executable `./test`