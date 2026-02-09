headers = $(wildcard inc/*.h)
src = $(wildcard src/*.c)
objs = matrix.o vector.o
blas_library = -lopenblas
CFLAGS ?= -Wall -fvectorize -fPIC -O2 -march=native
INSTALL_DIR ?= /usr/local

main:
	$(CC) -include $(headers) $(CCFLAGS) -c $(src)
	$(CC) -shared -o libmatrice.so $(objs)

test:
	@make main
	$(CC) $(objs) -O tests/test_allocations.c -o test_allocations
	$(CC) $(objs) -O tests/test_basics.c -o test_basics
	$(CC) $(objs) -O tests/test_ops.c -o test_ops
	$(CC) tests/test_shared.c -o test_shared libmatrice.so

test-cpp:
	@make main
	$(CXX) tests-cpp/test_bindings_cpp.cxx -o test_bindings_cpp libmatrice.so

bench:
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) -include $(headers) $(CCFLAGS) -c src/matrix.c -o matrix.o
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-float.c -o bench-float
	$(CC) matrix.o $(blas_library) -O benchmarks/bench-double.c -o bench-double

install:
	@echo Installing libmatrice.so
	cp libmatrice.so $(INSTALL_DIR)/lib
	@echo Installing header files
	mkdir -p $(INSTALL_DIR)/include/matrice
	cp inc/*.h $(INSTALL_DIR)/include/matrice

install-cpp:
	mkdir -p $(INSTALL_DIR)/include/matrice
	cp bindings/cpp/matrice.hxx $(INSTALL_DIR)/include/matrice

clean:
	rm -f *.o
	rm -f *.so
	rm -f inc/*.pch
	rm -f test_*

clean-bench:
	rm bench-float
	rm bench-double

