################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/algo_imr.cpp \
../src/main.cpp \
../src/matrix.cpp \
../src/vldb.cpp 

OBJS += \
./src/algo_imr.o \
./src/main.o \
./src/matrix.o \
./src/vldb.o 

CPP_DEPS += \
./src/algo_imr.d \
./src/main.d \
./src/matrix.d \
./src/vldb.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -I/opt/local/include -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


