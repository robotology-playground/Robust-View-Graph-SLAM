################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../SparseLib++/1.7/testing/tpre.cc \
../SparseLib++/1.7/testing/tsl.cc \
../SparseLib++/1.7/testing/tsp.cc \
../SparseLib++/1.7/testing/tspsm.cc 

OBJS += \
./SparseLib++/1.7/testing/tpre.o \
./SparseLib++/1.7/testing/tsl.o \
./SparseLib++/1.7/testing/tsp.o \
./SparseLib++/1.7/testing/tspsm.o 

CC_DEPS += \
./SparseLib++/1.7/testing/tpre.d \
./SparseLib++/1.7/testing/tsl.d \
./SparseLib++/1.7/testing/tsp.d \
./SparseLib++/1.7/testing/tspsm.d 


# Each subdirectory must supply rules for building sources it contributes
SparseLib++/1.7/testing/%.o: ../SparseLib++/1.7/testing/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


