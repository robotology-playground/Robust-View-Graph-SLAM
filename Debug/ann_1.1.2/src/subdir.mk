################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ann_1.1.2/src/ANN.cpp \
../ann_1.1.2/src/bd_fix_rad_search.cpp \
../ann_1.1.2/src/bd_pr_search.cpp \
../ann_1.1.2/src/bd_search.cpp \
../ann_1.1.2/src/bd_tree.cpp \
../ann_1.1.2/src/brute.cpp \
../ann_1.1.2/src/kd_dump.cpp \
../ann_1.1.2/src/kd_fix_rad_search.cpp \
../ann_1.1.2/src/kd_pr_search.cpp \
../ann_1.1.2/src/kd_search.cpp \
../ann_1.1.2/src/kd_split.cpp \
../ann_1.1.2/src/kd_tree.cpp \
../ann_1.1.2/src/kd_util.cpp \
../ann_1.1.2/src/perf.cpp 

OBJS += \
./ann_1.1.2/src/ANN.o \
./ann_1.1.2/src/bd_fix_rad_search.o \
./ann_1.1.2/src/bd_pr_search.o \
./ann_1.1.2/src/bd_search.o \
./ann_1.1.2/src/bd_tree.o \
./ann_1.1.2/src/brute.o \
./ann_1.1.2/src/kd_dump.o \
./ann_1.1.2/src/kd_fix_rad_search.o \
./ann_1.1.2/src/kd_pr_search.o \
./ann_1.1.2/src/kd_search.o \
./ann_1.1.2/src/kd_split.o \
./ann_1.1.2/src/kd_tree.o \
./ann_1.1.2/src/kd_util.o \
./ann_1.1.2/src/perf.o 

CPP_DEPS += \
./ann_1.1.2/src/ANN.d \
./ann_1.1.2/src/bd_fix_rad_search.d \
./ann_1.1.2/src/bd_pr_search.d \
./ann_1.1.2/src/bd_search.d \
./ann_1.1.2/src/bd_tree.d \
./ann_1.1.2/src/brute.d \
./ann_1.1.2/src/kd_dump.d \
./ann_1.1.2/src/kd_fix_rad_search.d \
./ann_1.1.2/src/kd_pr_search.d \
./ann_1.1.2/src/kd_search.d \
./ann_1.1.2/src/kd_split.d \
./ann_1.1.2/src/kd_tree.d \
./ann_1.1.2/src/kd_util.d \
./ann_1.1.2/src/perf.d 


# Each subdirectory must supply rules for building sources it contributes
ann_1.1.2/src/%.o: ../ann_1.1.2/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


