################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ann_1.1.2/sample/ann_sample.cpp 

OBJS += \
./ann_1.1.2/sample/ann_sample.o 

CPP_DEPS += \
./ann_1.1.2/sample/ann_sample.d 


# Each subdirectory must supply rules for building sources it contributes
ann_1.1.2/sample/%.o: ../ann_1.1.2/sample/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


