################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/5point.cpp \
../src/GraphOptimiser.cpp \
../src/PwgOptimiser.cpp \
../src/RecoverMoments.cpp \
../src/Rpoly.cpp \
../src/Tracker.cpp \
../src/VLFeat.cpp \
../src/main.cpp 

C_SRCS += \
../src/calibrated_fivepoint_helper.c 

OBJS += \
./src/5point.o \
./src/GraphOptimiser.o \
./src/PwgOptimiser.o \
./src/RecoverMoments.o \
./src/Rpoly.o \
./src/Tracker.o \
./src/VLFeat.o \
./src/calibrated_fivepoint_helper.o \
./src/main.o 

C_DEPS += \
./src/calibrated_fivepoint_helper.d 

CPP_DEPS += \
./src/5point.d \
./src/GraphOptimiser.d \
./src/PwgOptimiser.d \
./src/RecoverMoments.d \
./src/Rpoly.d \
./src/Tracker.d \
./src/VLFeat.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


