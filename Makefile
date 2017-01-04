all:
	g++ -o example -Wall -O2 example.C -lcrypto
	./example >example.out

clean:
	-rm example
