files=$(wildcard *.cpp)
main: $(patsubst %.cpp,  %.o, $(files))
	g++  $^ -O3 -o $@
%.o: %.cpp
	g++ -O3 -c $<
result: main *.h
	./main 
	 python3 *.py
	#Rscript graf*r
clear:
	rm -f *.o main
 
