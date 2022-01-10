#TODO: Get number of CUDA cores: https://stackoverflow.com/questions/32530604/how-can-i-get-number-of-cores-in-cuda-device

function gpuinfo()
    dev = device()
    println("                CUDA version: ", CUDA.CUDA_VERSION)
    println("        Multiprocessor count: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT))
    println("       Max threads per block: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK))
    println("     Max registers per block: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK))
    println(" Max shared memory per block: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK))
    println("          Max threads per MP: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR))
    println("        Max registers per MP: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_MULTIPROCESSOR))
    println("    Max shared memory per MP: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR))
    println("               L2 cache size: ", attribute(dev, CUDA.DEVICE_ATTRIBUTE_L2_CACHE_SIZE))
    println("                   warp size: ", CUDA.warpsize(dev))
    println("                Total memory: ", CUDA.total_memory())
    println("      Total memory on device: ", CUDA.totalmem(dev))
    println("          Compute capability: ", CUDA.capability(dev))
end

macro hascuda(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end