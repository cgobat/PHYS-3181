CFLAGS += -std=c99

generate_random_numbers: generate_random_numbers.o  mt19937-64.o
	gcc -o $@ $^ 
