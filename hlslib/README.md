_Do you use hlslib? Please consider [citing us](https://arxiv.org/abs/1910.04436), and let us know so we can feature your project in the list of examples._

## Quick introduction

hlslib is a collection of C++ headers, CMake files, and examples, aimed at improving the quality of life of HLS developers. The current repertoire primarily supports Vivado HLS, but some Intel FPGA OpenCL support is being added. An extended abstract describing the project is [available here](https://arxiv.org/abs/1910.04436).

This project is developed at the [Scalable Parallel Computing Lab](https://spcl.inf.ethz.ch/) (SPCL) at ETH Zurich (see our [github](https://github.com/spcl)).

#### How do I install it?

There are a few ways:
- Grab the headers and/or CMake files you need and stick them in your project.
- Install hlslib using the standard CMake installation procedure to a location of your choice.
- Clone this repository into your project as a git submodule and integrate it, with or without CMake.

#### How do I use it?

Just `#include` the header(s) you are interested in. You can see an example [here](https://github.com/spcl/gemm_hls).

When a Xilinx hlslib header is included, compilation must allow C++11 features, and the macro `HLSLIB_SYNTHESIS` must be set whenever HLS is run. Set `-cflags "-std=c++11 -DHLSLIB_SYNTHESIS"` in your synthesis script, and `--xp prop:kernel.<entry function>.kernel_flags="-std=c++11 -DHLSLIB_SYNTHESIS"` when building SDAccel kernels. See the included `xilinx_test/CMakeLists.txt` for reference. 

## Feature overview

We have Doxygen now! Simply run `make` to generate the docs.

A brief overview of hlslib features is given below.

#### CMake integration

For integrating the Xilinx or Intel HLS tools in your project, the `FindSDAccel.cmake` and `FindIntelFPGAOpenCL.cmake` are provided in the `cmake` subdirectory. The scripts will set all necessary variables required to build both host and device code.

Example `CMakeLists.txt`:
```cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} hlslib/cmake)
find_package(Vitis REQUIRED)

add_executable(MyHostExecutable src/MyHostExecutable.cpp)
include_directories(${Vitis_INCLUDE_DIRS})
target_link_libraries(MyHostExecutable ${Vitis_LIBRARIES})

add_custom_target(compile_kernel ${Vitis_COMPILER} -c -t hw src/MyKernel.cpp
                  --kernel MyKernel --platform xilinx_vcu1525_dynamic_5_1
                  -o MyKernel.xo)
add_custom_target(link_kernel ${Vitis_COMPILER} -l -t hw MyKernel.xo 
                  --kernel MyKernel --platform xilinx_vcu1525_dynamic_5_1
                  -o MyKernel.xclbin)
```

#### DataPack

The `hlslib::DataPack` class located in `hlslib/xilinx/DataPack.h` facilitates SIMD-style vectorization, and makes it easy to build wide data paths in your code.

Examples usage:
```cpp
hlslib::DataPack<float, 4> Foo(hlslib::DataPack<float, 4> &a,
                               hlslib::DataPack<float, 4> &b) {
  #pragma HLS INLINE     
  auto add = a + b; // Vector addition
  add[0] = add[1];  // Indexing for both reads and writes
  return 0.5 * add; // Element-wise multiplication by a scalar
}
```

#### Simulation

For kernels with multiple processing elements (PEs) executing in parallel, the `hlslib/xilinx/Simulation.h` adds some convenient macros to simulate this behavior, by wrapping each PE in a thread executed in parallel, all of which are joined when the program terminates.

Example usage:
```cpp
HLSLIB_DATAFLOW_INIT();
hlslib::Stream<Data_t> pipes[kStages + 1];
HLSLIB_DATAFLOW_FUNCTION(MemoryToStream, memory_in, pipes[0]);
for (int i = 0; i < kStages; ++i) {
  #pragma HLS UNROLL
  HLSLIB_DATAFLOW_FUNCTION(PE, pipes[i], pipes[i + 1]); // Launches new C++ thread
}
HLSLIB_DATAFLOW_FUNCTION(StreamToMemory, pipes[kStages], memory_out);
HLSLIB_DATAFLOW_FINALIZE(); // In simulation mode, joins threads created as dataflow functions.
```

When building programs using the simulation features, you must link against a thread library (e.g., pthreads).

#### Stream

While Vivado HLS provides the `hls::stream` class, it is somewhat lacking in features, in particular when simulating multiple processing elements. The `hlslib::Stream` class in `hlslib/xilinx/Stream.h` compiles to Vivado HLS streams, but provides a richer interface. hlslib streams are:
- thread-safe during simulation, allowing producer and consumer to be executed in parallel;
- bounded, simulating the finite capacity of hardware FIFOs, allowing easier detection of deadlocks in software; and
- self-contained, allowing the stream depth and implementation (e.g., using LUTRAM or BRAM) to be specified directly in the object, without excess pragmas.

Example usage:
```cpp
void Bar(hlslib::Stream<int> &a, hlslib::Stream<int> &b, int N) {
  for (int i = 0; i < N; ++i) {
    #pragma HLS PIPELINE II=1
    auto read = a.Pop(); // Queue-like interface
    b.Push(read + 1);
  }
}

void Foo(hlslib::Stream<int> &in_stream, // Specifying stream depth is optional
         hlslib::Stream<int> &out_stream, int N) {
  #pragma HLS DATAFLOW
  
  hlslib::Stream<int, 4> foo_pipe; // Implements a FIFO of depth 4
  
  // Dataflow functions running in parallel
  HLSLIB_DATAFLOW_INIT();
  HLSLIB_DATAFLOW_FUNCTION(Bar, in_stream, foo_pipe, N);
  HLSLIB_DATAFLOW_FUNCTION(Bar, foo_pipe, out_stream, N);
  HLSLIB_DATAFLOW_FINALIZE();
}
```

#### OpenCL host code

To greatly reduce the amount of boilerplate code required to create and launch OpenCL kernels, and to handle FPGA-specific configuration required by the vendors, hlslib provides a C++14 convenience interface in `hlslib/xilinx/SDAccel.h` and `hlslib/intel/OpenCL.h` for Vivado HLS and Intel FPGA OpenCL, respectively.

Example usage:
```cpp
using hlslib::ocl;
Context context;
std::vector<float> input_host(N, 5);   
std::vector<float> output_host(N, 5);
auto input_device = context.MakeBuffer<float, Access::read>(
    MemoryBank::bank0, input_host.cbegin(), input_end.cend());
auto output_device = context.MakeBuffer<float, Access::write>(MemoryBank::bank1, N);
auto program = context.MakeProgram("MyKernel.xclbin");
auto kernel = program.MakeKernel("MyKernel", input_device, output_device, N);
kernel.ExecuteTask();
output_device.CopyToHost(output_host.begin());
```

#### Other features

Various other features are provided, including:
* Classes to flatten loop nests and keep track of indices (`include/hlslib/xilinx/Flatten.h`), both for bounds known at runtime (`hlslib::Flatten`) and bounds known at compile-time (`hlslib::ConstFlatten`). Example usage can be found in `xilinx_test/kernels/Flatten.cpp`.
* Various compile-time functions commonly used when designing hardware, such as log2, in `include/hlslib/xilinx/Utility.h`.
* A template tcl-file that can be used with CMake or another templating engine to produce a high-level synthesis script.
* `xilinx_test/CMakeLists.txt` that builds a number of tests to verify hlslib functionality, doubling as a reference for how to integrate HLS projects with CMake using the provided module files .
* An example of how to use the Simulation and Stream headers, at `xilinx_test/kernels/MultiStageAdd.cpp`, both as a host-only simulation (`xilinx_test/test/TestMultiStageAdd.cpp`), and as a hardware kernel (`xilinx_test/host/RunMultiStageAdd.cpp`). 
* `include/hlslib/xilinx/Accumulate.h`, which includes a streaming implementation of accumulation, including for type/operator combinations with a non-zero latency (such as floating point addition). Example kernels of usage for both integer and floating point types are included as `xilinx_test/kernel/AccumulateInt.cpp` and `xilinx_test/kernel/AccumulateFloat.cpp`, respectively. 
* `include/hlslib/xilinx/Operators.h`, which includes some commonly used operators as functors to be plugged into templated functions such as `TreeReduce` and `Accumulate`.
* `include/hlslib/xilinx/Axi.h`, which implements the AXI Stream interface and the bus interfaces required by the DataMover IP, enabling the use of a command stream-based memory interface for HLS kernels if packaged as an RTL kernel where the DataMover IP is connected to the AXI interfaces.

Some of these headers depend on others. Please refer to the source code.

### Ubuntu packages

On Ubuntu, the following package might need to be installed to run hardware emulation:

```
sudo apt install libc6-dev-i386
```

## Projects using hlslib

- [Matrix multiplication code](https://github.com/spcl/gemm_hls) [1]: uses a wide range of hlslib features, including simulation, streams, vectors, CMake integration, and OpenCL wrapper code.
- [SMI](https://github.com/spcl/SMI) [2]: streaming message passing library for inter-FPGA communication in OpenCL. Uses hlslib for OpenCL host code.
- [HelmGemm](https://ieeexplore.ieee.org/document/8825124/) [3]: uses the simulation features of hlslib, and incorporates the matrix multiplication code above.
- [DaCe](https://arxiv.org/abs/1902.10345) [4]: a data-centric parallel programming framework targeting a wide range of architectures, including both Intel and Xilinx FPGAs. Uses hlslib for CMake integration, OpenCL host code, vectors, streams, and simulation.

_If you use hlslib in your project, please let us know, so we can add you to the list._

## Bugs and feature requests

Please use the issue tracker.

## Publication

If your project uses hlslib, please consider citing us:

**BibTeX:**
```
@article{hlslib,
  title={{hlslib}: Software Engineering for Hardware Design},
  author={de~Fine~Licht, Johannes and Hoefler, Torsten},
  journal={arXiv preprint arXiv:1910.04436},
  year={2019}
}
```

**Plain text:**
```
J. de Fine Licht and T. Hoefler, "hlslib: Software Engineering for Hardware Design", arXiv preprint arXiv:1910.04436 (2019).
```

## References

- [1] Johannes de Fine Licht, Grzegorz Kwasniewski, and Torsten Hoefler, _"Flexible Communication Avoiding Matrix Multiplication on FPGA with High-Level Synthesis"_, in Proceedings of 28th ACM/SIGDA International Symposium on Field-Programmable Gate Arrays (FPGA'20).
- [2] Tiziano De Matteis, Johannes de Fine Licht, Jakub Beránek, and Torsten Hoefler. _"Streaming Message Interface: High-Performance Distributed Memory Programming on Reconfigurable Hardware"_, in Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC'19).
- [3] Dionysios Diamantopoulos, and Christoph Hagleitner. _"HelmGemm: Managing GPUs and FPGAs for transprecision GEMM workloads in containerized environments."_, in Proceedings of the 2019 IEEE 30th International Conference on Application-specific Systems, Architectures and Processors (ASAP'19).
- [4] Tal Ben-Nun, Johannes de Fine Licht, Alexandros Nikolaos Ziogas, Timo Schneider, and Torsten Hoefler. _"Stateful Dataflow Multigraphs: A Data-Centric Model for High-Performance Parallel Programs"_, in Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC'19).
