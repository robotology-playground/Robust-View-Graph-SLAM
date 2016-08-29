################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../mex/mex_compute_gate_graph.cpp \
../mex/mex_compute_gate_inverse_depth_Mviews.cpp \
../mex/mex_constraints_addition_inverse_depth_Mviews.cpp \
../mex/mex_generate_constraints_info_Mviews.cpp \
../mex/mex_observation_model_inverse_depth.cpp \
../mex/mex_observation_model_inverse_depth_Mviews.cpp \
../mex/mex_observation_model_jacobian_inverse_depth.cpp \
../mex/mex_observation_model_jacobian_inverse_depth_Mviews.cpp \
../mex/mex_optimise_constraints_Mviews.cpp \
../mex/mex_recover_moments.cpp \
../mex/mex_update_info_matrix_Mviews.cpp 

OBJS += \
./mex/mex_compute_gate_graph.o \
./mex/mex_compute_gate_inverse_depth_Mviews.o \
./mex/mex_constraints_addition_inverse_depth_Mviews.o \
./mex/mex_generate_constraints_info_Mviews.o \
./mex/mex_observation_model_inverse_depth.o \
./mex/mex_observation_model_inverse_depth_Mviews.o \
./mex/mex_observation_model_jacobian_inverse_depth.o \
./mex/mex_observation_model_jacobian_inverse_depth_Mviews.o \
./mex/mex_optimise_constraints_Mviews.o \
./mex/mex_recover_moments.o \
./mex/mex_update_info_matrix_Mviews.o 

CPP_DEPS += \
./mex/mex_compute_gate_graph.d \
./mex/mex_compute_gate_inverse_depth_Mviews.d \
./mex/mex_constraints_addition_inverse_depth_Mviews.d \
./mex/mex_generate_constraints_info_Mviews.d \
./mex/mex_observation_model_inverse_depth.d \
./mex/mex_observation_model_inverse_depth_Mviews.d \
./mex/mex_observation_model_jacobian_inverse_depth.d \
./mex/mex_observation_model_jacobian_inverse_depth_Mviews.d \
./mex/mex_optimise_constraints_Mviews.d \
./mex/mex_recover_moments.d \
./mex/mex_update_info_matrix_Mviews.d 


# Each subdirectory must supply rules for building sources it contributes
mex/%.o: ../mex/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


