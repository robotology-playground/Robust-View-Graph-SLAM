################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../mex/scratch/mex_SCC.cpp \
../mex/scratch/mex_compute_gate_inverse_depth.cpp \
../mex/scratch/mex_formulate_system.cpp \
../mex/scratch/mex_kdeSample.cpp \
../mex/scratch/mex_least_squares.cpp \
../mex/scratch/mex_mfas.cpp \
../mex/scratch/mex_observation_model.cpp \
../mex/scratch/mex_observation_model_jacobian.cpp \
../mex/scratch/mex_recover_moments.cpp \
../mex/scratch/mex_recover_second_moment.cpp \
../mex/scratch/mex_rotation_averaging.cpp \
../mex/scratch/mex_solve_translation.cpp \
../mex/scratch/mex_triangulate.cpp 

C_SRCS += \
../mex/scratch/phonebook.c 

OBJS += \
./mex/scratch/mex_SCC.o \
./mex/scratch/mex_compute_gate_inverse_depth.o \
./mex/scratch/mex_formulate_system.o \
./mex/scratch/mex_kdeSample.o \
./mex/scratch/mex_least_squares.o \
./mex/scratch/mex_mfas.o \
./mex/scratch/mex_observation_model.o \
./mex/scratch/mex_observation_model_jacobian.o \
./mex/scratch/mex_recover_moments.o \
./mex/scratch/mex_recover_second_moment.o \
./mex/scratch/mex_rotation_averaging.o \
./mex/scratch/mex_solve_translation.o \
./mex/scratch/mex_triangulate.o \
./mex/scratch/phonebook.o 

C_DEPS += \
./mex/scratch/phonebook.d 

CPP_DEPS += \
./mex/scratch/mex_SCC.d \
./mex/scratch/mex_compute_gate_inverse_depth.d \
./mex/scratch/mex_formulate_system.d \
./mex/scratch/mex_kdeSample.d \
./mex/scratch/mex_least_squares.d \
./mex/scratch/mex_mfas.d \
./mex/scratch/mex_observation_model.d \
./mex/scratch/mex_observation_model_jacobian.d \
./mex/scratch/mex_recover_moments.d \
./mex/scratch/mex_recover_second_moment.d \
./mex/scratch/mex_rotation_averaging.d \
./mex/scratch/mex_solve_translation.d \
./mex/scratch/mex_triangulate.d 


# Each subdirectory must supply rules for building sources it contributes
mex/scratch/%.o: ../mex/scratch/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

mex/scratch/%.o: ../mex/scratch/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


