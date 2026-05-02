
test_eq: src/eq.cpp src/sim.cpp src/address_eval.cpp src/test_eq.cpp include/eq.hpp src/thal.cpp include/thal.h
	g++ -Wall src/eq.cpp src/sim.cpp src/address_eval.cpp src/test_eq.cpp src/thal.cpp -Iinclude -O3 -lm -lpthread -mavx -o test_eq

clean:
	rm test_eq
