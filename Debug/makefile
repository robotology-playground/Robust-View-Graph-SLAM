################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/imageviewer/subdir.mk
-include src/subdir.mk
-include mex/scratch/subdir.mk
-include mex/subdir.mk
-include build/CMakeFiles/vgSLAM.dir/src/subdir.mk
-include build/CMakeFiles/3.6.1/CompilerIdCXX/subdir.mk
-include build/CMakeFiles/3.6.1/CompilerIdC/subdir.mk
-include build/CMakeFiles/3.2.2/CompilerIdCXX/subdir.mk
-include build/CMakeFiles/3.2.2/CompilerIdC/subdir.mk
-include build/CMakeFiles/subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.managedbuilder.core/subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: Robust-View-Graph-SLAM

# Tool invocations
Robust-View-Graph-SLAM: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++  -o "Robust-View-Graph-SLAM" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) Robust-View-Graph-SLAM
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
