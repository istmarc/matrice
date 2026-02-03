CC = clang

headers = $(wildcard inc/*.h)
src = $(wildcard src/*.c)

main:
	$(CC) -include $(headers) -c $(src)

test:
	$(CC) -include $(headers) -c src/matrix.c -o matrix.o
	$(CC) matrix.o -O tests/test_allocations.c -o test_allocations


