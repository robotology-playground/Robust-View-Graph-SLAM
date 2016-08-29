################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ann_1.1.2/ann2fig/ann2fig.cpp 

OBJS += \
./ann_1.1.2/ann2fig/ann2fig.o 

CPP_DEPS += \
./ann_1.1.2/ann2fig/ann2fig.d 


# Each subdirectory must supply rules for building sources it contributes
ann_1.1.2/ann2fig/%.o: ../ann_1.1.2/ann2fig/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


