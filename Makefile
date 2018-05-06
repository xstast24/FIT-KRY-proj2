default: clean
	g++ main.cpp -o kry -lgmp # lgmp is GMP lib
	./kry -g 8

clean:
	rm -f *.o
	rm -f kry
