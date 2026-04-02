#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <CL/cl2.hpp>
#include <iostream>

int main() {
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    if (all_platforms.size()==0) {
        std::cout << "No platforms found. Check OpenCL installation\n";
        exit(1);
    }
    cl::Platform default_platform=all_platforms[0];
    std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if (all_devices.size()=0) {
        std::cout << " No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Device default_device=all_devices[0];
    std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";

    cl::Context context({default_device});
    cl::Program::Sources sources;

    std::string kernal_code=
        " void kernal simple_add(global const int* A, global const int* B, global int* C){ "
        "   C[get_global_id(0)]=A[get_global_id(0)]+B[get_global_id(0)];    "
        " } ";
    sources.push_back({kernal_code.c_str(), kernal_code.length()});

    cl::Program program(context, sources);
    if (program.build({default_device}) != CL_SUCCESS) {
        std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    cl::Buffer buffer_A(context, CL_MEM_READ_WRITE, sizeof(int)*10);
    cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, sizeof(int)*10);
    cl::Buffer buffer_C(context, CL_MEM_READ_WRITE, sizeof(int)*10);

    int A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};

    cl::CommandQueue queue(context, default_device);

    queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(int)*10, A);
    queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(int)*10, B);


    cl::KernelFunctor simple_add(cl::Kernel(program,"simple_add"),queue,cl::NullRange,cl::NDRange(10),cl::NullRange);
    simple_add(buffer_A, buffer_B, buffer_C);

    int C[10];
    queue.enqueueReadBuffer(buffer_C, CL_TRUE, 0, sizeof(int)*10, C);

    std::cout << "results: " << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << C[i] << " ";
    }

    return 0;
}
