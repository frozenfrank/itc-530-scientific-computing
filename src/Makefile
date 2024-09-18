.PHONY: clean
wave_solve: main.cpp wave_orthotope.hpp wave_orthotope.cpp
	g++ -std=c++20 -o $@ $< wave_orthotope.cpp

serial_debug: main_serial.cpp wave_orthotope.hpp wave_orthotope.cpp binary_io.hpp
	g++ -std=c++20 -g -o $@ $< wave_orthotope.cpp ../src/binary_io.hpp
clean:
	$(RM) wave_solve
