################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../SparseLib++/1.7/mv/src/mvblasc.cc \
../SparseLib++/1.7/mv/src/mvblasd.cc \
../SparseLib++/1.7/mv/src/mvblasf.cc \
../SparseLib++/1.7/mv/src/mvblasi.cc \
../SparseLib++/1.7/mv/src/mvblast.cc \
../SparseLib++/1.7/mv/src/mvmc.cc \
../SparseLib++/1.7/mv/src/mvmd.cc \
../SparseLib++/1.7/mv/src/mvmf.cc \
../SparseLib++/1.7/mv/src/mvmi.cc \
../SparseLib++/1.7/mv/src/mvmt.cc \
../SparseLib++/1.7/mv/src/mvvc.cc \
../SparseLib++/1.7/mv/src/mvvcio.cc \
../SparseLib++/1.7/mv/src/mvvd.cc \
../SparseLib++/1.7/mv/src/mvvdio.cc \
../SparseLib++/1.7/mv/src/mvvf.cc \
../SparseLib++/1.7/mv/src/mvvi.cc \
../SparseLib++/1.7/mv/src/mvvt.cc 

OBJS += \
./SparseLib++/1.7/mv/src/mvblasc.o \
./SparseLib++/1.7/mv/src/mvblasd.o \
./SparseLib++/1.7/mv/src/mvblasf.o \
./SparseLib++/1.7/mv/src/mvblasi.o \
./SparseLib++/1.7/mv/src/mvblast.o \
./SparseLib++/1.7/mv/src/mvmc.o \
./SparseLib++/1.7/mv/src/mvmd.o \
./SparseLib++/1.7/mv/src/mvmf.o \
./SparseLib++/1.7/mv/src/mvmi.o \
./SparseLib++/1.7/mv/src/mvmt.o \
./SparseLib++/1.7/mv/src/mvvc.o \
./SparseLib++/1.7/mv/src/mvvcio.o \
./SparseLib++/1.7/mv/src/mvvd.o \
./SparseLib++/1.7/mv/src/mvvdio.o \
./SparseLib++/1.7/mv/src/mvvf.o \
./SparseLib++/1.7/mv/src/mvvi.o \
./SparseLib++/1.7/mv/src/mvvt.o 

CC_DEPS += \
./SparseLib++/1.7/mv/src/mvblasc.d \
./SparseLib++/1.7/mv/src/mvblasd.d \
./SparseLib++/1.7/mv/src/mvblasf.d \
./SparseLib++/1.7/mv/src/mvblasi.d \
./SparseLib++/1.7/mv/src/mvblast.d \
./SparseLib++/1.7/mv/src/mvmc.d \
./SparseLib++/1.7/mv/src/mvmd.d \
./SparseLib++/1.7/mv/src/mvmf.d \
./SparseLib++/1.7/mv/src/mvmi.d \
./SparseLib++/1.7/mv/src/mvmt.d \
./SparseLib++/1.7/mv/src/mvvc.d \
./SparseLib++/1.7/mv/src/mvvcio.d \
./SparseLib++/1.7/mv/src/mvvd.d \
./SparseLib++/1.7/mv/src/mvvdio.d \
./SparseLib++/1.7/mv/src/mvvf.d \
./SparseLib++/1.7/mv/src/mvvi.d \
./SparseLib++/1.7/mv/src/mvvt.d 


# Each subdirectory must supply rules for building sources it contributes
SparseLib++/1.7/mv/src/%.o: ../SparseLib++/1.7/mv/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


