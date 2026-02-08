headers = $(wildcard inc/*.h)
src = $(wildcard src/*.c)
blas_library = -lopenblas
CCFLAGS = -Wall -fvectorize -O2 -march=native

main:
	$(CC) -include $(headers) -c $(src)

test:
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) matrix.o -O tests/test_allocations.c -o test_allocations
	$(CC) matrix.o -O tests/test_basics.c -o test_basics
	$(CC) matrix.o -O tests/test_ops.c -o test_ops

bench:
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-float.c -o bench-float
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-double.c -o bench-double



