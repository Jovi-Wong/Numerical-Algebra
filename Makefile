source = $(wildcard *.cpp)
#object = $(patsubst %.cpp, %.o, $(source))

main: $(object)
	g++ -o $@ $(source)

debug: $(source)
	g++ -o main $(source) -g

clean:
	rm *.o main 
	rm -rf main.dSYM
