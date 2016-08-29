################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../SparseLib++/1.7/src/iohb.c \
../SparseLib++/1.7/src/iotext.c 

CC_SRCS += \
../SparseLib++/1.7/src/compcol_double.cc \
../SparseLib++/1.7/src/comprow_double.cc \
../SparseLib++/1.7/src/coord_double.cc \
../SparseLib++/1.7/src/diagpre.cc \
../SparseLib++/1.7/src/diagpre_double.cc \
../SparseLib++/1.7/src/icpre.cc \
../SparseLib++/1.7/src/icpre_double.cc \
../SparseLib++/1.7/src/ilupre.cc \
../SparseLib++/1.7/src/ilupre_1.5.cc \
../SparseLib++/1.7/src/ilupre_double.cc \
../SparseLib++/1.7/src/ilupre_double_1.5.cc \
../SparseLib++/1.7/src/iohb_double.cc \
../SparseLib++/1.7/src/iotext.cc \
../SparseLib++/1.7/src/iotext_double.cc \
../SparseLib++/1.7/src/qsort_double.cc \
../SparseLib++/1.7/src/qsort_int.cc \
../SparseLib++/1.7/src/qsort_type.cc 

OBJS += \
./SparseLib++/1.7/src/compcol_double.o \
./SparseLib++/1.7/src/comprow_double.o \
./SparseLib++/1.7/src/coord_double.o \
./SparseLib++/1.7/src/diagpre.o \
./SparseLib++/1.7/src/diagpre_double.o \
./SparseLib++/1.7/src/icpre.o \
./SparseLib++/1.7/src/icpre_double.o \
./SparseLib++/1.7/src/ilupre.o \
./SparseLib++/1.7/src/ilupre_1.5.o \
./SparseLib++/1.7/src/ilupre_double.o \
./SparseLib++/1.7/src/ilupre_double_1.5.o \
./SparseLib++/1.7/src/iohb.o \
./SparseLib++/1.7/src/iohb_double.o \
./SparseLib++/1.7/src/iotext.o \
./SparseLib++/1.7/src/iotext_double.o \
./SparseLib++/1.7/src/qsort_double.o \
./SparseLib++/1.7/src/qsort_int.o \
./SparseLib++/1.7/src/qsort_type.o 

C_DEPS += \
./SparseLib++/1.7/src/iohb.d \
./SparseLib++/1.7/src/iotext.d 

CC_DEPS += \
./SparseLib++/1.7/src/compcol_double.d \
./SparseLib++/1.7/src/comprow_double.d \
./SparseLib++/1.7/src/coord_double.d \
./SparseLib++/1.7/src/diagpre.d \
./SparseLib++/1.7/src/diagpre_double.d \
./SparseLib++/1.7/src/icpre.d \
./SparseLib++/1.7/src/icpre_double.d \
./SparseLib++/1.7/src/ilupre.d \
./SparseLib++/1.7/src/ilupre_1.5.d \
./SparseLib++/1.7/src/ilupre_double.d \
./SparseLib++/1.7/src/ilupre_double_1.5.d \
./SparseLib++/1.7/src/iohb_double.d \
./SparseLib++/1.7/src/iotext.d \
./SparseLib++/1.7/src/iotext_double.d \
./SparseLib++/1.7/src/qsort_double.d \
./SparseLib++/1.7/src/qsort_int.d \
./SparseLib++/1.7/src/qsort_type.d 


# Each subdirectory must supply rules for building sources it contributes
SparseLib++/1.7/src/%.o: ../SparseLib++/1.7/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

SparseLib++/1.7/src/%.o: ../SparseLib++/1.7/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


