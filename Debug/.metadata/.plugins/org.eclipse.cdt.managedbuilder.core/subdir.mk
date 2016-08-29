################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).C \
../.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec.C 

OBJS += \
./.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).o \
./.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec.o 

C_UPPER_DEPS += \
./.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).d \
./.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec.d 


# Each subdirectory must supply rules for building sources it contributes
.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).o: ../.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).C
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF".metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec (Case Conflict).d" -MT".metadata/.plugins/org.eclipse.cdt.managedbuilder.core/spec\ (Case\ Conflict).d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/%.o: ../.metadata/.plugins/org.eclipse.cdt.managedbuilder.core/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


