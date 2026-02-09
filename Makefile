headers = $(wildcard inc/*.h)
src = $(wildcard src/*.c)
objs = matrix.o vector.o
blas_library = -lopenblas
CCFLAGS = -Wall -fvectorize -fPIC -O2 -march=native

main:
	$(CC) -include $(headers) $(CCFLAGS) -c $(src)
	$(CC) -shared -o matrice.so $(objs)

test:
	$(CC) -include $(headers) $(CCFLAGS) -c $(src)
	$(CC) -shared -o matrice.so $(objs)
	$(CC) $(objs) -O tests/test_allocations.c -o test_allocations
	$(CC) $(objs) -O tests/test_basics.c -o test_basics
	$(CC) $(objs) -O tests/test_ops.c -o test_ops
	$(CC) tests/test_shared.c -o test_shared matrice.so

bench:
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-float.c -o bench-float
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-double.c -o bench-double

clean:
	rm *.o
	rm *.so
	rm inc/*.pch

clean-tests:
	rm test_*

clean-bench:
	rm bench-float
	rm bench-double

