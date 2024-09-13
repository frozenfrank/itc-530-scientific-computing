.PHONY: clean
wave_solve: main.cpp wave_orthotope.hpp wave_orthotope.cpp
	g++ -std=c++20 -o $@ $< wave_orthotope.cpp

clean:
	$(RM) wave_solve
