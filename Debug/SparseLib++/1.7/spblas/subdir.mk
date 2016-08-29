################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../SparseLib++/1.7/spblas/spmm.cc \
../SparseLib++/1.7/spblas/spsm.cc 

OBJS += \
./SparseLib++/1.7/spblas/spmm.o \
./SparseLib++/1.7/spblas/spsm.o 

CC_DEPS += \
./SparseLib++/1.7/spblas/spmm.d \
./SparseLib++/1.7/spblas/spsm.d 


# Each subdirectory must supply rules for building sources it contributes
SparseLib++/1.7/spblas/%.o: ../SparseLib++/1.7/spblas/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


