default: clean
	g++ main.cpp -o kry -lgmp # lgmp is GMP lib

clean:
	rm -f *.o
	rm -f kry
