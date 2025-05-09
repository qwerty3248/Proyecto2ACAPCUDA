NVCC = nvcc
CC = gcc
SRC_DIR = src
BIN_DIR = bin

# Archivos fuente
CUDA_SRC = $(SRC_DIR)/FourierCUDA.cu
SEQ_SRC = $(SRC_DIR)/FourierS.c

# Ejecutables
CUDA_EXE = $(BIN_DIR)/FourierCUDA
SEQ_EXE = $(BIN_DIR)/FourierS


all: $(CUDA_EXE) $(SEQ_EXE)

$(CUDA_EXE): $(CUDA_SRC)
	$(NVCC) -std=c++11 $< -o $@ -lcudart -lm

$(SEQ_EXE): $(SEQ_SRC)
	$(CC) -std=c11 $< -o $@ -lm

clean:
	rm -rf $(BIN_DIR)/*
